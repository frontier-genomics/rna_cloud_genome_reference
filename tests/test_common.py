import tempfile
import pytest
import os

from rnacloud_genome_reference.grc_fixes.common import download_file

class TestCommon:
    def test_successful_download_file(self):
        try:
            temp_dir = tempfile.mkdtemp()
            temp_file = tempfile.mkstemp()
            download_file('https://raw.githubusercontent.com/github/.github/refs/heads/main/profile/README.md', temp_dir, temp_file[1])
        except Exception as e:
            pytest.fail(f"download_file raised an exception: {e}")
        finally:
            os.remove(temp_file[1])
            os.rmdir(temp_dir)

    def test_unsuccessful_download_file(self):
        with pytest.raises(Exception):
            download_file('https://raw.com/github/.github/refs/heads/main/profile/README.md', temp_dir, temp_file[1])
    
