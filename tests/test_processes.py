"""Cleanup failures must not replace an assembly/annotation failure."""
import signal
import subprocess
from unittest.mock import Mock, patch

import pytest

from basic_assembler import processes


def child():
    process = Mock(spec=subprocess.Popen)
    process.pid = 12345
    process.stdin = Mock()
    process.stdout = Mock()
    process.stderr = Mock()
    return process


@pytest.mark.parametrize('error', [PermissionError('denied'), ProcessLookupError('gone')])
def test_group_signal_failure_falls_back_to_direct_child(error):
    process = child()
    with patch.object(processes.os, 'name', 'posix'), patch.object(processes.os, 'killpg', side_effect=error) as killpg:
        processes.stop_process_tree(process)
    killpg.assert_called_once_with(process.pid, signal.SIGKILL)
    process.kill.assert_called_once_with()


def test_group_is_signalled_even_after_wrapper_exits():
    process = child()
    process.returncode = 0
    process.poll.return_value = 0
    with patch.object(processes.os, 'name', 'posix'), patch.object(processes.os, 'killpg') as killpg:
        processes.stop_process_tree(process)
    killpg.assert_called_once_with(process.pid, signal.SIGKILL)
    process.kill.assert_not_called()


@pytest.mark.parametrize('wait_error', [subprocess.TimeoutExpired('aligner', 1), OSError('wait failed')])
def test_denied_signals_use_bounded_wait_and_close_pipes(wait_error):
    process = child()
    process.kill.side_effect = PermissionError('PID denied')
    process.wait.side_effect = wait_error
    with patch.object(processes.os, 'name', 'posix'), patch.object(processes.os, 'killpg', side_effect=PermissionError('group denied')):
        processes.cleanup_process(process)
    process.wait.assert_called_once_with(timeout=1.0)
    for stream in (process.stdin, process.stdout, process.stderr):
        stream.close.assert_called_once_with()
    process.communicate.assert_not_called()


def test_close_failures_preserve_original_error_and_other_pipes_close():
    process = child()
    process.stdin.close.side_effect = OSError('close failed')
    process.stdout.close.side_effect = ValueError('detached stream')
    original = RuntimeError('original assembly failure')
    with patch.object(processes.os, 'name', 'posix'), patch.object(processes.os, 'killpg'):
        with pytest.raises(RuntimeError) as caught:
            try:
                raise original
            finally:
                processes.cleanup_process(process)
    assert caught.value is original
    process.stderr.close.assert_called_once_with()


def test_non_posix_stops_direct_child_and_tolerates_uncaptured_pipes():
    process = child()
    process.stdin = process.stdout = process.stderr = None
    with patch.object(processes.os, 'name', 'nt'), patch.object(processes.os, 'killpg') as killpg:
        processes.cleanup_process(process)
    killpg.assert_not_called()
    process.kill.assert_called_once_with()
    process.wait.assert_called_once_with(timeout=1.0)
