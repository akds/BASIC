"""Best-effort subprocess cleanup that preserves the caller's original error."""
from __future__ import annotations

import os
import signal
import subprocess


_CLEANUP_TIMEOUT = 1.0


def stop_process_tree(process: subprocess.Popen) -> None:
    """Signal descendants even after the wrapper exits; fall back to its PID."""
    if os.name == 'posix':
        try:
            os.killpg(process.pid, signal.SIGKILL)
            return
        except OSError:
            pass
    try:
        process.kill()
    except OSError:
        pass


def cleanup_process(process: subprocess.Popen) -> None:
    """Stop and briefly reap a failed subprocess, then close captured pipes.

    Do not call this from a watchdog while another thread reads the pipes;
    watchdogs should use ``stop_process_tree`` instead. Never drain inherited
    pipes with communicate(), since surviving descendants may keep them open.
    """
    stop_process_tree(process)
    try:
        process.wait(timeout=_CLEANUP_TIMEOUT)
    except (OSError, subprocess.TimeoutExpired):
        pass
    for stream in (process.stdin, process.stdout, process.stderr):
        if stream is not None:
            try:
                stream.close()
            except (OSError, ValueError):
                pass
