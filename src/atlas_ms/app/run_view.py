"""
Running the pipeline from the app: Run / Stop, progress and log.

The app never runs the tools itself (see CLAUDE.md): Run starts the same
Snakemake command as ``atlas-ms run`` in a subprocess, and its output is
read line by line without blocking the interface. Snakemake reports
"N of M steps (P%) done" after each step, which drives the progress bar.
"""

import asyncio
import re
from collections import deque
from typing import Callable

import panel as pn

from atlas_ms.project import Project
from atlas_ms.runner import snakemake_command

PROGRESS = re.compile(r"(\d+) of (\d+) steps \((\d+)%\) done")
LOG_LINES = 40  # lines of Snakemake output shown under the buttons


class RunPanel:
    """Sidebar controls that run the workflow on the open project."""

    def __init__(
        self,
        project: Project,
        before_run: Callable[[], bool],
        after_run: Callable[[], None],
        command: Callable[..., list[str]] = snakemake_command,
    ):
        """
        ``before_run`` saves the project and returns False to cancel (invalid
        input); ``after_run`` reloads the results after a successful run;
        ``command`` builds the Snakemake command (replaced in tests).
        """
        self.project = project
        self.before_run, self.after_run, self.command = before_run, after_run, command
        self.process = None
        self.lines = deque(maxlen=LOG_LINES)

        self.cores = pn.widgets.IntInput(label="CPU cores", value=4, start=1, width=120)
        self.run_button = pn.widgets.Button(label="Run", color="primary", width=100)
        self.stop_button = pn.widgets.Button(label="Stop", color="danger", width=100, disabled=True)
        self.progress = pn.indicators.Progress(value=0, max=100, sizing_mode="stretch_width")
        self.status = pn.pane.Markdown("Ready.")
        self.log = pn.pane.HTML("", sizing_mode="stretch_width")
        self.run_button.on_click(self.run)
        self.stop_button.on_click(self.stop)
        self.layout = pn.Column(
            self.cores, pn.Row(self.run_button, self.stop_button), self.progress, self.status, self.log,
            sizing_mode="stretch_width",
        )

    def _show(self, line: str) -> None:
        self.lines.append(line)
        text = "\n".join(self.lines).replace("&", "&amp;").replace("<", "&lt;")
        self.log.object = f"<pre style='font-size:11px; white-space:pre-wrap'>{text}</pre>"
        match = PROGRESS.search(line)
        if match:
            self.progress.value = int(match.group(3))
            self.status.object = f"Step {match.group(1)} of {match.group(2)}"

    async def run(self, event=None) -> int | None:
        """Save, run the workflow, show its output; returns Snakemake's exit code."""
        if self.process is not None or not self.before_run():
            return None
        self.lines.clear()
        self.progress.value, self.status.object = 0, "Running..."
        self.run_button.disabled, self.stop_button.disabled = True, False
        try:
            self.process = await asyncio.create_subprocess_exec(
                *self.command(self.project.root, cores=self.cores.value),
                stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.STDOUT,
            )
            while line := await self.process.stdout.readline():
                self._show(line.decode(errors="replace").rstrip())
            code = await self.process.wait()
        finally:
            self.process = None
            self.run_button.disabled, self.stop_button.disabled = False, True
        if code == 0:
            self.progress.value, self.status.object = 100, "Finished."
            self.after_run()
        else:
            self.status.object = f"Failed (exit code {code}): see the log below and the files in logs/."
        return code

    def stop(self, event=None) -> None:
        """Stop Snakemake; interrupted steps are redone on the next run."""
        if self.process is not None:
            self.process.terminate()
            self.status.object = "Stopping..."
