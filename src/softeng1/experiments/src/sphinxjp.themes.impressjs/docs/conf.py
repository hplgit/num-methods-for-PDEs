# -*- coding: utf-8 -*-
#
# -- General configuration -------------------------------------

source_suffix = '.rst'
master_doc = 'index'

project = u'Sphinx theme for dynamic html presentation style'
copyright = u'2012, Sphinx-users.jp'

version = '0.1.2'

# -- Options for HTML output -----------------------------------

extensions = ['sphinxjp.themecore',
              'sphinxcontrib.blockdiag',
              'sphinxcontrib.seqdiag']
html_theme = 'impressjs'
html_use_index = False
