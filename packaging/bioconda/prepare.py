#!/usr/bin/env python3
"""Generate a local Bioconda recipe from an actual BASIC source distribution."""
from __future__ import annotations

import argparse
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
from string import Template
import tarfile
import tempfile
from urllib.parse import urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener


HERE = Path(__file__).resolve().parent


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _inspect_sdist(path: Path, license_file: str) -> dict:
    license_path = PurePosixPath(license_file)
    if license_path.is_absolute() or '..' in license_path.parts or not license_file:
        raise ValueError('license_file must be a relative path within the sdist')
    with tarfile.open(path, 'r:gz') as archive:
        members = archive.getmembers()
        names = set()
        for member in members:
            parts = PurePosixPath(member.name)
            if (not parts.parts or parts.is_absolute() or '..' in parts.parts
                    or not (member.isfile() or member.isdir())):
                raise ValueError(f'Unsafe source archive member: {member.name}')
            if member.name in names:
                raise ValueError(f'Duplicate source archive member: {member.name}')
            names.add(member.name)
        metadata = [member for member in members
                    if len(PurePosixPath(member.name).parts) == 2 and member.name.endswith('/PKG-INFO')]
        if len(metadata) != 1 or not metadata[0].isfile() or metadata[0].size > 2 * 1024 * 1024:
            raise ValueError('Expected one root PKG-INFO file of at most 2 MiB')
        root = PurePosixPath(metadata[0].name).parts[0]
        if any(PurePosixPath(member.name).parts[0] != root for member in members):
            raise ValueError('Source archive must contain one top-level directory')
        for required in ('pyproject.toml', 'basic_assembler/__init__.py', license_file):
            member = archive.getmember(f'{root}/{required}') if f'{root}/{required}' in names else None
            if member is None or not member.isfile() or member.size == 0:
                raise ValueError(f'Source archive is missing a nonempty {required}')
        handle = archive.extractfile(metadata[0])
        assert handle is not None
        with handle:
            info = BytesParser().parsebytes(handle.read())
    name, version = info.get('Name', ''), info.get('Version', '')
    if re.sub(r'[-_.]+', '-', name).lower() != 'basic-immune-assembler':
        raise ValueError(f'Unexpected distribution name: {name!r}')
    if not re.fullmatch(r'[0-9][A-Za-z0-9.+]*', version):
        raise ValueError(f'Invalid or missing source version: {version!r}')
    return {'distribution': name, 'version': version}


class _HTTPSRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urlsplit(newurl).scheme != 'https':
            raise ValueError('Source URL redirected outside HTTPS')
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _verify_public_source(url: str, expected_sha256: str, expected_size: int) -> None:
    parsed = urlsplit(url)
    if (parsed.scheme != 'https' or not parsed.hostname or parsed.username is not None
            or parsed.password is not None or parsed.fragment):
        raise ValueError('source_url must be a public HTTPS URL without credentials or a fragment')
    parsed.port  # Reject invalid ports before making a request.
    digest, size = hashlib.sha256(), 0
    request = Request(url, headers={'User-Agent': 'BASIC-Bioconda-preparation/1'})
    with build_opener(_HTTPSRedirect()).open(request, timeout=30) as response:
        while chunk := response.read(1024 * 1024):
            size += len(chunk)
            if size > expected_size:
                raise ValueError('Published source is larger than the local sdist')
            digest.update(chunk)
    if size != expected_size or digest.hexdigest() != expected_sha256:
        raise ValueError('Published source does not match the local sdist SHA256 and size')


def prepare(sdist: Path, output: Path, *, source_url: str | None = None,
            license_expression: str = 'MIT',
            license_file: str = 'LICENSE') -> dict:
    sdist, output = sdist.expanduser().resolve(), output.expanduser().absolute()
    if output.exists():
        raise FileExistsError(f'Recipe output already exists: {output}; choose a new directory')
    if not license_expression or len(license_expression) > 200 or any(ord(c) < 32 for c in license_expression):
        raise ValueError('license_expression must be a nonempty, single-line expression')
    info = _inspect_sdist(sdist, license_file)
    info.update(sdist=str(sdist), sha256=_sha256(sdist), bytes=sdist.stat().st_size,
                source_url=source_url or sdist.as_uri(),
                source_mode='verified_https' if source_url else 'local_file',
                license_expression=license_expression, license_file=license_file)
    if source_url:
        _verify_public_source(source_url, info['sha256'], info['bytes'])
    template = Template((HERE / 'meta.yaml.in').read_text(encoding='utf-8'))
    recipe = template.substitute({key: json.dumps(info[key]) for key in
                                  ('version', 'source_url', 'sha256', 'license_expression', 'license_file')})
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.basic-recipe-', dir=output.parent) as temporary:
        stage = Path(temporary) / 'recipe'
        stage.mkdir()
        (stage / 'meta.yaml').write_text(recipe, encoding='utf-8')
        shutil.copy2(HERE / 'run_test.py', stage / 'run_test.py')
        (stage / 'preparation.json').write_text(json.dumps(info, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        stage.rename(output)
    return info


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sdist', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--source-url', help='Existing HTTPS URL; downloaded bytes must match the local sdist')
    parser.add_argument('--license-expression', default='MIT')
    parser.add_argument('--license-file', default='LICENSE')
    args = parser.parse_args(argv)
    try:
        info = prepare(**vars(args))
    except (OSError, ValueError, tarfile.TarError) as exc:
        parser.exit(1, f'Bioconda preparation failed: {exc}\n')
    print(json.dumps(info, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
