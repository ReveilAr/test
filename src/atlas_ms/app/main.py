"""
The ATLAS-MS app: one page with a sidebar (project, run) and tabs.

    atlas-ms app [project folder] [--port 5006]

* Sidebar: open an existing project or create one from a folder of raw
  files; Run / Stop the pipeline, with progress and log.
* Setup tab: sample metadata and parameters (``setup_view``).
* Network tab: molecular network, feature list, spectrum and chromatograms
  (``network_view``), once the pipeline has produced results.
* Annotation tab: every annotation candidate of the selected feature and
  its evidence (``annotation_view``).
* Quit (bottom of the sidebar): stops a running pipeline and the app
  itself, which frees its port (an app left running keeps its port busy).

One project is open at a time. Opening another one rebuilds the page.
"""

import socket
from pathlib import Path

import panel as pn

from atlas_ms.app.annotation_view import AnnotationTab
from atlas_ms.app.data import load_results
from atlas_ms.preprocessing.msdata import RunReader
from atlas_ms.app.network_view import NetworkTab
from atlas_ms.app.run_view import RunPanel
from atlas_ms.app.setup_view import SetupTab
from atlas_ms.app.uploads import drop_area
from atlas_ms.config import load_presets
from atlas_ms.project import RAW_FORMATS, Project
from atlas_ms.runner import snakemake_command

pn.extension("tabulator")

ACCENT = "#c51b8a"  # magenta


class AtlasApp:
    """The page, and the project currently open in it."""

    def __init__(self, project_dir: str | None = None, command=snakemake_command, on_quit=None):
        """
        ``command`` builds the Snakemake command (replaced in tests);
        ``on_quit`` stops the server (given by ``serve``; no Quit button without it).
        """
        self.command = command
        self.on_quit = on_quit
        self.project = None
        self.run_panel = None
        self.content = pn.Column(sizing_mode="stretch_width")  # the tabs of the open project
        self.run_area = pn.Column(sizing_mode="stretch_width")  # the Run panel of the open project

        # ---- open / create a project ----
        self.path = pn.widgets.TextInput(label="Project folder", value=project_dir or "", sizing_mode="stretch_width")
        open_button = pn.widgets.Button(label="Open", color="primary", width=100)
        open_button.on_click(lambda event: self.open(self.path.value))
        self.raw_dir = pn.widgets.TextInput(label="Raw files folder (.raw / .mzML)", sizing_mode="stretch_width")
        # Or drop the files: they are copied into <project folder>/raw/ (see uploads.py).
        self.raw_drop = drop_area("", [".raw", ".mzml"], "...or drop .raw / .mzML files here (copied into the project)")
        self.raw_drop.param.watch(self._raw_rejected, "rejected")
        self.path.param.watch(self._project_folder_changed, "value")
        self._project_folder_changed()
        self.new_instrument = pn.widgets.Select(label="Instrument", options=list(load_presets("instruments")))
        self.new_adducts = pn.widgets.Select(label="Adducts", options=list(load_presets("adducts")))
        create_button = pn.widgets.Button(label="Create project", width=140)
        create_button.on_click(lambda event: self.create())

        self.quit_button = pn.widgets.Button(label="Quit ATLAS-MS", color="danger", width=160,
                                             visible=on_quit is not None)
        self.quit_button.on_click(lambda event: self.quit())

        self.sidebar = [
            self.path, open_button,
            pn.Card(self.raw_dir, self.raw_drop, self.new_instrument, self.new_adducts, create_button,
                    title="New project (in the folder above)", collapsed=True, sizing_mode="stretch_width"),
            pn.layout.Divider(),
            self.run_area,
            pn.layout.Divider(),
            self.quit_button,
        ]
        if project_dir:
            self.open(project_dir)
        else:
            self.content.objects = [pn.pane.Markdown("Open a project folder, or create one, in the sidebar.")]

    def open(self, folder: str) -> None:
        """Open a project: build its Setup tab, Run panel and (if there are results) Network tab."""
        project = Project(folder)
        if not project.config_path.exists():
            self.content.objects = [pn.pane.Alert(f"No ATLAS-MS project in {project.root}", alert_type="danger")]
            return
        self.project = project
        self.setup = SetupTab(project)
        self.run_panel = RunPanel(project, before_run=self.setup.save, after_run=self.refresh, command=self.command)
        self.run_area.objects = [self.run_panel.layout]
        self.refresh()

    def _project_folder_changed(self, event=None) -> None:
        """Dropped raw files go into <project folder>/raw/ (nowhere until a folder is typed)."""
        folder = self.path.value.strip()
        self.raw_drop.folder = str(Path(folder).expanduser() / "raw") if folder else ""

    def _raw_rejected(self, event) -> None:
        reason = "type the project folder first" if not self.raw_drop.folder else "not a .raw or .mzML file"
        self.content.objects = [pn.pane.Alert(f"Not added ({reason}): {', '.join(event.new[len(event.old):])}",
                                              alert_type="warning")]

    def create(self) -> None:
        """
        Create a project in the "Project folder" from every raw file of the raw
        files folder and every dropped file.
        """
        files = []
        if self.raw_dir.value.strip():
            raw_dir = Path(self.raw_dir.value).expanduser()
            files = sorted(f for f in raw_dir.iterdir() if f.suffix.lower() in RAW_FORMATS) if raw_dir.is_dir() else []
        files += [Path(f) for f in self.raw_drop.saved if Path(f) not in files]
        if not files:
            self.content.objects = [pn.pane.Alert("No .raw or .mzML file: give a raw files folder or drop files",
                                                  alert_type="danger")]
            return
        try:
            Project.create(self.path.value, files, instrument=self.new_instrument.value, adducts=self.new_adducts.value)
        except (FileExistsError, ValueError) as error:
            self.content.objects = [pn.pane.Alert(str(error), alert_type="danger")]
            return
        self.open(self.path.value)

    def refresh(self) -> None:
        """(Re)build the tabs from the files on disk, e.g. after a run."""
        tabs = [("Setup", self.setup.layout)]
        results = load_results(self.project)
        if results is not None:
            reader = RunReader(self.project.root, list(self.project.load_samples().index))
            tolerance = self.setup.config.scoring.fragment_tolerance_da
            self.network = NetworkTab(results, reader, fragment_tolerance=tolerance)
            self.annotation = AnnotationTab(results, self.network)
            tabs.append(("Network", self.network.view()))
            tabs.append(("Annotation", self.annotation.view()))
        else:
            tabs.append(("Network", pn.pane.Markdown("No results yet: press **Run** in the sidebar.")))
        self.tabs = pn.Tabs(*tabs, sizing_mode="stretch_width", dynamic=True)
        if results is not None:
            self.tabs.active = 1
        self.content.objects = [self.tabs]

    def quit(self) -> None:
        """
        First click: ask for a confirmation. Second click: stop a running
        pipeline, say goodbye and stop the server (``on_quit``).
        """
        if self.quit_button.label != "Click again to quit":
            self.quit_button.label = "Click again to quit"
            return
        if self.run_panel is not None:
            self.run_panel.stop()
        self.content.objects = [pn.pane.Alert("ATLAS-MS has stopped: you can close this tab. "
                                              "Start it again with `atlas-ms app`.", alert_type="info")]
        self.on_quit()

    def page(self) -> pn.template.BaseTemplate:
        return pn.template.FastListTemplate(
            title="ATLAS-MS", sidebar=self.sidebar, main=[self.content], sidebar_width=340,
            accent=ACCENT,  # header bar and primary buttons
        )


def free_port(port: int, attempts: int = 20) -> int:
    """
    The first port, from ``port`` on, that nothing is listening on. A port
    is often still taken by an app left running in another terminal.
    """
    for candidate in range(port, port + attempts):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
            try:
                probe.bind(("", candidate))  # all interfaces, as the app server does
            except OSError:
                continue
            return candidate
    raise OSError(f"Ports {port}-{port + attempts - 1} are all in use: choose another one with --port")


def stop_later(server, delay_s: float = 1.0) -> None:
    """
    Stop the web server after ``delay_s`` (time for the page to show its
    goodbye): the port is freed and ``serve`` returns, ending the process.
    """
    def stop():
        server.stop()
        server.io_loop.stop()
    server.io_loop.call_later(delay_s, stop)


def serve(project_dir: str | None = None, port: int = 5006, show: bool = True) -> None:
    """Start the app server; each browser tab gets its own page. Returns when Quit is clicked."""
    chosen = free_port(port)
    if chosen != port:
        print(f"Port {port} is already in use (another ATLAS-MS app still running?): using port {chosen}.")
    server = pn.serve(lambda: AtlasApp(project_dir, on_quit=lambda: stop_later(server)).page(),
                      port=chosen, show=show, title="ATLAS-MS", start=False)
    server.start()
    try:
        server.io_loop.start()  # runs until stop_later (Quit) or Ctrl+C
    except KeyboardInterrupt:
        pass
    print("ATLAS-MS stopped.")
