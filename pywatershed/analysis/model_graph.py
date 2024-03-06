import pathlib as pl
import tempfile

from ..base.conservative_process import ConservativeProcess
from ..base.model import Model
from ..utils import import_optional_dependency


class ModelGraph:
    def __init__(
        self,
        model: Model,
        show_params: bool = False,
        process_colors: dict = None,
        node_penwidth: int = 2,
        default_edge_color: str = "black",
        from_file_edge_color: str = None,
        node_spacing: float = 2.75,
        hide_variables: bool = True,
    ):
        self.pydot = import_optional_dependency("pydot")
        self.graph = None

        self.model = model
        self.show_params = show_params
        self.process_colors = process_colors
        self.node_penwidth = node_penwidth
        self.default_edge_color = default_edge_color
        if not from_file_edge_color:
            self.from_file_edge_color = default_edge_color
        else:
            self.from_file_edge_color = from_file_edge_color
        self.node_spacing = node_spacing
        self.hide_variables = hide_variables

        return

    def build_graph(self):
        # Build the process nodes in the graph
        self._current_pos = self.node_spacing
        self.process_nodes = {}
        for process in self.model.process_order:
            self.process_nodes[process] = self._process_node(
                process,
                self.model.processes[process],
                show_params=self.show_params,
            )

        # Solve the connections
        self.files = []
        self.connections = []
        for process in self.model.process_order:
            frm_already = []
            for var, frm in self.model.process_input_from[process].items():
                var_con = f":{var}"

                if self.hide_variables:
                    var_con = ""
                    if frm in frm_already:
                        # print(f"skipping: {frm}")
                        continue
                    else:
                        frm_already += [frm]
                        # print(frm_already)

                if not isinstance(frm, pl.Path):
                    color = self.default_edge_color
                    if self.process_colors:
                        color = self.process_colors[frm]
                    self.connections += [
                        (
                            f"{frm}{var_con}",
                            f"{process}{var_con}",
                            color,
                        )
                    ]

                else:
                    file_name = frm.name
                    self.files += [file_name]
                    self.connections += [
                        (
                            f"Files:{file_name.split('.')[0]}",
                            f"{process}{var_con}",
                            self.from_file_edge_color,
                        )
                    ]

        # Build the file node, reset the position
        self._current_pos = 0
        self.file_node = self._file_node(self.files)

        # build the graph
        self.graph = self.pydot.Dot(
            graph_type="digraph",
            layout="neato",
            splines="polyline",
        )

        self.graph.add_node(self.file_node)

        for process in self.model.process_order:
            self.graph.add_node(self.process_nodes[process])

        for con in self.connections:
            self.graph.add_edge(self.pydot.Edge(con[0], con[1], color=con[2]))

        return

    def SVG(self, verbose: bool = False, dpi=45):
        """Display an SVG in jupyter notebook (via tempfile)."""

        ipdisplay = import_optional_dependency("IPython.display")

        tmp_file = pl.Path(tempfile.NamedTemporaryFile().name)
        if self.graph is None:
            self.build_graph()
        self.graph.write_svg(tmp_file, prog=["dot", f"-Gdpi={dpi}"])
        if verbose:
            print(f"Displaying SVG written to temp file: {tmp_file}")

        ipdisplay.display(ipdisplay.SVG(tmp_file))
        return

    def _process_node(self, process_name, process, show_params: bool = False):
        inputs = process.get_inputs()
        variables = process.get_variables()
        params = process.get_parameters()
        cls = process.__class__.__name__

        label = (
            f'<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="1">\n'
            f'    <TR><TD COLSPAN="6">"{process_name}"</TD></TR>\n'
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

        if isinstance(process, ConservativeProcess):
            mass_budget_vars = [
                var
                for comp, vars in process.get_mass_budget_terms().items()
                for var in vars
            ]
        else:
            mass_budget_vars = []

        for varset_name, varset in show_categories.items():
            n_vars = len(varset)
            if not n_vars:
                continue

            varset_col_span = 2
            if self.hide_variables:
                varset = []
                n_vars = 0
                varset_col_span = 6

            label += f"    <TR>\n"
            label += (
                f'        <TD COLSPAN="{varset_col_span}" ROWSPAN="{n_vars+1}"'
            )
            label += f'         BGCOLOR="{category_colors[varset_name]}">{varset_name}</TD>\n'
            label += f"    </TR>\n"

            for ii, vv in enumerate(sorted(varset)):
                border_color_str = ""
                if vv in mass_budget_vars:
                    border_color_str = 'border="1" COLOR="BLUE"'
                label += f"    <TR>\n"
                label += f'        <TD COLSPAN="4" BGCOLOR="{category_colors[varset_name]}" {border_color_str} PORT="{vv}"><FONT POINT-SIZE="9.0">{vv}</FONT></TD>\n'
                label += f"    </TR>\n"

        label += f"</TABLE>>\n"
        label = label

        color_str = ""
        if self.process_colors:
            color_str = f'"{self.process_colors[process_name]}"'

        node = self.pydot.Node(
            process_name,
            label=label,
            pos=f'"{self._current_pos},0!"',
            shape="box",
            color=color_str,
            penwidth=f'"{self.node_penwidth}"',
        )
        self._current_pos += self.node_spacing

        return node

    def _file_node(self, files):
        files = list(set(files))
        label = (
            f'<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="1">\n'
            f'    <TR><TD COLSPAN="1">Files</TD></TR>\n'
        )

        for file in files:
            label += f"    <TR>\n"
            label += f'        <TD COLSPAN="1"  BGCOLOR="gray50" PORT="{file.split(".")[0]}"><FONT POINT-SIZE="9.0">{file}</FONT></TD>\n'
            label += f"    </TR>\n"

        label += f"</TABLE>>\n"
        label = label
        node = self.pydot.Node(
            "Files",
            label=label,
            pos=f'"{self._current_pos},0!"',
            shape="note",
            color=f'"{self.from_file_edge_color}"',
            penwidth=f'"{self.node_penwidth}"',
        )
        self._current_pos += self.node_spacing
        return node
