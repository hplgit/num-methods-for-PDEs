"""
This module contains a few small sphinx extensions.
They are mainly used to help with the generation
of BPS's own documentation, but some other projects
use them as well, so they are kept here.
"""
#=============================================================================
# imports
#=============================================================================
# core
import re
import os.path
import sys
# local
__all__ = [
    # constants
    "__version__",
    "std_exts",
    "all_exts",

    # public utility functions
    "get_theme_dir",
    "get_version",

    # internal helpers
    "is_cloud_theme",
    "u",
]

#=============================================================================
# constants
#=============================================================================
__version__ = "1.6"

# names of standard cloud extensions
# used by most cloud themes
std_exts = [
    'cloud_sptheme.ext.autodoc_sections',
    'cloud_sptheme.ext.index_styling',
    'cloud_sptheme.ext.relbar_toc',
    'cloud_sptheme.ext.table_styling',
]

# names of all cloud extensions
all_exts = std_exts + [
    'cloud_sptheme.ext.issue_tracker',
    'cloud_sptheme.ext.escaped_samp_literals',
]

#=============================================================================
# public helpers
#=============================================================================
_root_dir = os.path.abspath(os.path.dirname(__file__))

def get_theme_dir():
    """Returns path to directory containing this package's Sphinx themes.

    This is designed to be used when setting the ``html_theme_path``
    option within Sphinx's ``conf.py`` file.

    .. seealso:: The :ref:`Cloud Sphinx Theme <cloud-theme-usage>` for a usage example.
    """
    return os.path.join(_root_dir, "themes")

def get_version(release):
    """Derive short version string from longer 'release' string.

    This is quick helper which takes a project's ``release`` string,
    and generates the shortened ``version`` string required by ``conf.py``.
    Usage example for ``conf.py``::

        import cloud_sptheme as csp

        ...

        # The version info for the project you're documenting
        from myapp import __version__ as release
        version = csp.get_version(release)
    """
    return re.match("(\d+\.\d+)", release).group(1)

#=============================================================================
# misc internal helpers
#=============================================================================
def is_cloud_theme(name):
    """[hack] internal helper to check if theme accepts cloud theme options"""
    return os.path.isfile(os.path.join(get_theme_dir(), name, "theme.conf"))

#=============================================================================
# internal py2/3 compat helpers
#=============================================================================
PY2 = sys.version_info < (3,0)
PY3 = not PY2
if PY2:
    def u(s):
        return s.decode("unicode_escape")
    def ru(s):
        return s.decode("ascii")
else:
    def u(s):
        return s
    ru = u

#=============================================================================
# eof
#=============================================================================
