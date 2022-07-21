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
        node_penwidth: int = 5,
    ):
        if not has_pydot:
            warnings.warn("pydot not available")

        self.model = model
        self.show_params = show_params
        self.component_colors = component_colors
        self.node_penwidth = node_penwidth
        self.graph = None

        return

    def build_graph(self):
        # Build the component nodes in the graph
        self.component_nodes = {}
        for component in self.model.component_order:
            self.component_nodes[component] = self._component_node(
                self.model.components[component], show_params=self.show_params
            )

        self.graph = pydot.Dot(
            graph_type="digraph", rankdir="LR", colorscheme="accent8"
        )

        for component in self.model.component_order:
            self.graph.add_node(self.component_nodes[component])

        # Downward connections
        for component in self.model.component_order:
            for var, frm in self.model.component_input_from[component].items():
                if not isinstance(frm, pl.Path):
                    self.graph.add_edge(
                        pydot.Edge(
                            f"{frm}:{var}",
                            f"{component}:{var}",
                            color=self.component_colors[frm],
                        )
                    )

        # "Upwards" connections across time.
        # graph.add_edge(pydot.Edge("PRMSSnow:pkwater_ante", "PRMSCanopy:pkwater_ante", color="blue"))

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

        # unclear correct syntax to set node directly
        # node = pydot.Node(label)
        # node.set_shape('plain')
        # go through a digraph instead
        color_str = ""
        if self.component_colors:
            color_str = f"""color="{self.component_colors[cls]}" """

        dot_string = (
            f"""digraph G {{ {cls} ["""
            f"""shape=component """
            f"""{color_str} """
            f"""penwidth={self.node_penwidth} """
            f"""label={label} """
            f"""]}}"""
        )
        node = pydot.graph_from_dot_data(dot_string)[0].get_node(cls)[0]

        return node
