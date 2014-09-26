# -*- coding: utf-8 -*-
"""
    sphinxjp.themes.imprssjs.directives
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Impress.js Style specific HTML tag builder.

    :copyright: Copyright (c) 2012, by Sphinx-users.jp, See AUTHORS.
    :license: MIT, see LICENSE for details.
"""

from docutils.parsers.rst import directives
from docutils.parsers.rst.roles import set_classes
from docutils import nodes
from sphinx.util.compat import Directive

__docformat__ = 'reStrructuredText'


class impressjs(nodes.General, nodes.Element): pass


class Impressjs(Directive):
    """
    A impressjs entry, control impressjs slide effects, actions and styles.
    """
    has_content = True
    required_arguments = 1
    optional_arguments = 1
    final_argument_whitespace = True

    option_spec = {
        'data-x': int,
        'data-y': int,
        'data-z': int,
        'data-rotate-x': int,
        'data-rotate-y': int,
        'data-rotate': int,
        'data-scale-x': int,
        'data-scale-y': int,
        'data-scale': int,
        'class': directives.class_option ,
    }

    node_class = impressjs

    def run(self):
        """ build impressjs node """
        set_classes(self.options)
        self.assert_has_content()
        text = '\n'.join(self.content)
        node = self.node_class(text, **self.options)
        self.add_name(node)
        self.state.nested_parse(self.content, self.content_offset, node)

        if self.arguments[0]:
            node['ids'] += [self.arguments[0]]
        if 'class' in self.options:
            node['classes'].append(options['class'])

        return [node]


def visit_impressjs(self, node):
    """ build div start tag for impres.js """
    atts = {'class': 'step'}
        
    if 'data-x' in node:
        atts['data-x'] = node['data-x']
    if 'data-y' in node:
        atts['data-y'] = node['data-y']
    if 'data-z' in node:
        atts['data-z'] = node['data-z']
    if 'data-rotate-x' in node:
        atts['data-rotate-x'] = node['data-rotate-x']
    if 'data-rotate-y' in node:
        atts['data-rotate-y'] = node['data-rotate-y']
    if 'data-rotate' in node:
        atts['data-rotate'] = node['data-rotate']
    if 'data-scale-x' in node:
        atts['data-scale-x'] = node['data-scale-x']
    if 'data-scale-y' in node:
        atts['data-scale-y'] = node['data-scale-y']
    if 'data-scale' in node:
        atts['data-scale'] = node['data-scale']

    self.body.append(self.starttag(node, 'div', **atts))
    self.set_first_last(node)


def depart_impressjs(self, node=None):
    """ build div end tag """
    self.body.append('</div>\n')


def setup(app):
    app.add_node(impressjs,
                 html=(visit_impressjs, depart_impressjs))
    app.add_directive('impressjs', Impressjs)
