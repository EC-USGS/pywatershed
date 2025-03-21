# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
from datetime import datetime

sys.path.insert(0, os.path.abspath("../"))
sys.path.insert(0, os.path.abspath("../pywatershed"))


import sphinx_autosummary_accessors  # noqa

import pywatershed  # noqa

# -- Project information -----------------------------------------------------

project = "pywatershed"
copyright = datetime.now().strftime("%Y")
author = "USGS Developers and Community"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# For sphinx-autodoc-typehints compat with sphinx.ext.napoleon see
# https://github.com/tox-dev/sphinx-autodoc-typehints/issues/15
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    # "nbsphinx",
    "sphinx_autosummary_accessors",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    # "sphinx_copybutton",
    # "sphinx.ext.mathjax",
    # "numpydoc",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
]

extlinks = {
    "issue": ("https://github.com/EC-USGS/pywatershed/issues/%s", "GH%s"),
    "pull": ("https://github.com/EC-USGS/pywatershed/pull/%s", "PR%s"),
}


autosummary_generate = True
# autodoc_typehints = "none"
# autosummary_imported_members = True

autodoc_default_options = {
    "members": True,
    # "imported-members": True,
    "inherited-members": True,
    "undoc-members": True,
    "private-members": False,  #
    # "special-members": "",
    "exclude-members": "__init__",
}


# Napoleon configurations
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_preprocess_types = True


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# The master toctree document.
master_doc = "index"

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"

html_title = "pywatershed"

html_context = {
    "github_user": "EC-USGS",
    "github_repo": "pywatershed",
    "github_version": "main",
    "doc_path": "doc",
}

html_logo = "_static/logo.png"
html_favicon = "_static/logo.ico"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["style.css"]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = dict(
    # analytics_id=''  this is configured in rtfd.io
    # canonical_url="",
    repository_url="https://github.com/EC-USGS/pywatershed",
    repository_branch="main",
    path_to_docs="doc",
    use_edit_page_button=True,
    use_repository_button=True,
    use_issues_button=True,
    # home_page_in_toc=True,
    # show_navbar_depth=1,
    # show_toc_level=1,
)

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "numba": ("https://numba.readthedocs.io/en/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "cftime": ("https://unidata.github.io/cftime", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
}
