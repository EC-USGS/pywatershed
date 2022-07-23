import pathlib as pl
import tempfile
import warnings

from ..base.model import Model

try:
    import pydot

    has_pydot = True
except ModuleNotFoundError:
    has_pydot = False

try:
    from IPython.display import SVG, display

    has_ipython = True
except ModuleNotFoundError:
    has_ipython = False


class ModelGraph:
    def __init__(
        self,
        model: Model,
        show_params: bool = False,
        component_colors: dict = None,
        node_penwidth: int = 2,
        default_edge_color: str = "gray50",
        from_file_edge_color: str = None,
        node_spacing=2.75,
    ):
        if not has_pydot:
            warnings.warn("pydot not available")

        self.model = model
        self.show_params = show_params
        self.component_colors = component_colors
        self.node_penwidth = node_penwidth
        self.default_edge_color = default_edge_color
        if not from_file_edge_color:
            self.from_file_edge_color = default_edge_color
        else:
            self.from_file_edge_color = from_file_edge_color
        self.node_spacing = node_spacing

        self.graph = None

        return

    def build_graph(self):
        # Build the component nodes in the graph
        self._current_pos = self.node_spacing
        self.component_nodes = {}
        for component in self.model.component_order:
            self.component_nodes[component] = self._component_node(
                self.model.components[component], show_params=self.show_params
            )

        # Solve the connections
        files = []
        connections = []
        for component in self.model.component_order:
            for var, frm in self.model.component_input_from[component].items():
                if not isinstance(frm, pl.Path):
                    color = self.default_edge_color
                    if self.component_colors:
                        color = self.component_colors[frm]
                    connections += [
                        (f"{frm}:{var}", f"{component}:{var}", color)
                    ]
                else:
                    file_name = frm.name
                    files += [file_name]
                    connections += [
                        (
                            f"Files:{file_name.split('.')[0]}",
                            f"{component}:{var}",
                            self.from_file_edge_color,
                        )
                    ]

        # Build the file node, reset the position
        self._current_pos = 0
        self.file_node = self._file_node(files)

        # build the graph
        self.graph = pydot.Dot(
            graph_type="digraph",
            layout="neato",
            splines="polyline",
        )

        self.graph.add_node(self.file_node)

        for component in self.model.component_order:
            self.graph.add_node(self.component_nodes[component])

        for con in connections:
            self.graph.add_edge(pydot.Edge(con[0], con[1], color=con[2]))

        return

    def SVG(self, verbose: bool = False):
        """Display an SVG in jupyter notebook (via tempfile)."""

        if not has_ipython:
            warnings.warn("IPython is not available")
        tmp_file = pl.Path(tempfile.NamedTemporaryFile().name)
        if self.graph is None:
            self.build_graph()
        self.graph.write_svg(tmp_file)
        if verbose:
            print(f"Displaying SVG written to temp file: {tmp_file}")

        display(SVG(tmp_file))
        return

    def _component_node(self, component, show_params: bool = False):
        inputs = component.get_inputs()
        variables = component.get_variables()
        params = component.get_parameters()
        cls = component.__class__.__name__

        label = (
            f'<<TABLE BORDER="0" CELLBORDER=".5" CELLSPACING="0" CELLPADDING="1">\n'
            f'    <TR><TD COLSPAN="6">{cls}</TD></TR>\n'
        )

        category_colors = {
            "inputs": "lightblue",
            "params": "orange",
            "variables": "lightgreen",
        }

        show_categories = {
            "inputs": inputs,
            "params": params,
            "variables": variables,
        }
        if not show_params:
            _ = show_categories.pop("params")

        for varset_name, varset in show_categories.items():

            n_vars = len(varset)
            if not n_vars:
                continue

            label += f"    <TR>\n"
            label += f'        <TD COLSPAN="2" ROWSPAN="{n_vars+1}"'
            label += f'         BGCOLOR="{category_colors[varset_name]}">{varset_name}</TD>\n'
            label += f"    </TR>\n"
            for ii, vv in enumerate(sorted(varset)):
                label += f"    <TR>\n"
                label += f'        <TD COLSPAN="4"  BGCOLOR="{category_colors[varset_name]}" PORT="{vv}"><FONT POINT-SIZE="9.0">{vv}</FONT></TD>\n'
                label += f"    </TR>\n"

        label += f"</TABLE>>\n"
        label = label

        color_str = ""
        if self.component_colors:
            color_str = f'"{self.component_colors[cls]}"'

        node = pydot.Node(
            cls,
            label=label,
            pos=f'"{self._current_pos},0!"',
            shape="component",
            color=color_str,
            penwidth=f'"{self.node_penwidth}"',
        )
        self._current_pos += self.node_spacing

        return node

    def _file_node(self, files):
        files = list(set(files))
        label = (
            f'<<TABLE BORDER="0" CELLBORDER=".5" CELLSPACING="0" CELLPADDING="1">\n'
            f'    <TR><TD COLSPAN="1">Files</TD></TR>\n'
        )

        for file in files:
            label += f"    <TR>\n"
            label += f'        <TD COLSPAN="1"  BGCOLOR="gray50" PORT="{file.split(".")[0]}"><FONT POINT-SIZE="9.0">{file}</FONT></TD>\n'
            label += f"    </TR>\n"

        label += f"</TABLE>>\n"
        label = label
        node = pydot.Node(
            "Files",
            label=label,
            pos=f'"{self._current_pos},0!"',
            shape="note",
            color=f'"{self.from_file_edge_color}"',
            penwidth=f'"{self.node_penwidth}"',
        )
        self._current_pos += self.node_spacing
        return node
