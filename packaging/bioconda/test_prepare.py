"""Standard-library tests for the recipe generator (no conda installation needed)."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location('basic_conda_prepare', Path(__file__).with_name('prepare.py'))
recipe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recipe)


class RecipeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def sdist(self, *, name='basic-immune-assembler', version='2.0.0', omit=None, extra=None):
        path = self.root / 'source release.tar.gz'
        files = {'PKG-INFO': f'Metadata-Version: 2.4\nName: {name}\nVersion: {version}\n',
                 'pyproject.toml': '[build-system]\n', 'basic_assembler/__init__.py': '__version__="2.0.0"\n',
                 'LICENSE': 'MIT License\nCopyright BASIC contributors\n'}
        if omit:
            del files[omit]
        with tarfile.open(path, 'w:gz') as archive:
            for name, contents in files.items():
                data = contents.encode()
                member = tarfile.TarInfo('source/' + name)
                member.size = len(data)
                archive.addfile(member, io.BytesIO(data))
            if extra:
                archive.addfile(extra)
        return path

    def test_generates_hashed_local_recipe_from_metadata(self):
        source = self.sdist(version='2.0.1rc1')
        output = self.root / 'recipe'
        info = recipe.prepare(source, output)
        self.assertEqual(info['version'], '2.0.1rc1')
        self.assertEqual(info['sha256'], hashlib.sha256(source.read_bytes()).hexdigest())
        self.assertEqual(info['source_url'], source.resolve().as_uri())
        self.assertEqual(info['license_expression'], 'MIT')
        self.assertEqual(info['license_file'], 'LICENSE')
        self.assertEqual(json.loads((output / 'preparation.json').read_text()), info)
        text = (output / 'meta.yaml').read_text()
        for expected in ('name: basic', 'python >=3.10', '- bowtie2', 'basic-references =', 'basic-skill =', info['sha256']):
            self.assertIn(expected, text)
        self.assertNotIn('github.com/akds/BASIC/archive/v2', text)
        self.assertTrue((output / 'run_test.py').is_file())

    def test_invalid_distribution_and_missing_license_do_not_create_recipe(self):
        for kwargs in ({'name': 'other-project'}, {'version': 'not-a-version'}, {'omit': 'LICENSE'}):
            with self.subTest(kwargs=kwargs):
                source = self.sdist(**kwargs)
                output = self.root / 'recipe'
                with self.assertRaises(ValueError):
                    recipe.prepare(source, output)
                self.assertFalse(output.exists())

    def test_existing_output_is_untouched(self):
        source = self.sdist()
        output = self.root / 'recipe'
        output.mkdir()
        marker = output / 'keep'
        marker.write_text('original')
        with self.assertRaises(FileExistsError):
            recipe.prepare(source, output)
        self.assertEqual(marker.read_text(), 'original')

    def test_unsafe_members_are_rejected_without_extraction(self):
        for name, kind in (('../outside', tarfile.REGTYPE), ('source/link', tarfile.SYMTYPE), ('.', tarfile.DIRTYPE)):
            with self.subTest(name=name):
                member = tarfile.TarInfo(name)
                member.type = kind
                member.linkname = '/outside'
                source = self.sdist(extra=member)
                with self.assertRaises(ValueError):
                    recipe.prepare(source, self.root / 'recipe')
        self.assertFalse((self.root / 'recipe').exists())

    def test_public_source_is_verified_before_recipe_creation(self):
        source = self.sdist()
        payload = source.read_bytes()
        with patch.object(recipe, 'build_opener') as factory:
            factory.return_value.open.return_value = io.BytesIO(payload)
            info = recipe.prepare(source, self.root / 'recipe', source_url='https://files.example.org/source.tar.gz')
        self.assertEqual(info['source_mode'], 'verified_https')
        self.assertEqual(info['source_url'], 'https://files.example.org/source.tar.gz')

    def test_public_source_mismatch_or_extra_bytes_leave_no_recipe(self):
        source = self.sdist()
        data = source.read_bytes()
        for payload in (b'wrong', bytes([data[0] ^ 1]) + data[1:], data + b'extra'):
            with self.subTest(length=len(payload)), patch.object(recipe, 'build_opener') as factory:
                factory.return_value.open.return_value = io.BytesIO(payload)
                with self.assertRaises(ValueError):
                    recipe.prepare(source, self.root / 'recipe', source_url='https://files.example.org/source.tar.gz')
            self.assertFalse((self.root / 'recipe').exists())

    def test_non_https_and_credentialed_urls_do_not_open_network(self):
        source = self.sdist()
        for url in ('http://example.org/source', 'https://user:password@example.org/source', 'https://example.org:bad/source'):
            with self.subTest(url=url), patch.object(recipe, 'build_opener') as factory:
                with self.assertRaises(ValueError):
                    recipe.prepare(source, self.root / 'recipe', source_url=url)
                factory.assert_not_called()

    def test_redirect_to_http_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'outside HTTPS'):
            recipe._HTTPSRedirect().redirect_request(None, None, 302, 'Found', {}, 'http://example.org/source')

    def test_cli_uses_archive_metadata(self):
        source = self.sdist()
        with patch('sys.stdout', new_callable=io.StringIO) as stdout:
            self.assertEqual(recipe.main(['--sdist', str(source), '--output', str(self.root / 'recipe')]), 0)
        self.assertEqual(json.loads(stdout.getvalue())['version'], '2.0.0')


if __name__ == '__main__':
    unittest.main()
