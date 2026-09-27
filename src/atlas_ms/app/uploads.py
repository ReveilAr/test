"""
Drag and drop of files into the app.

A browser never tells a web page where a dropped file is on the disk, only
what it contains. Dropping a file therefore *uploads* it: the app receives
it in pieces and writes a copy into a folder (the project's ``raw/`` or
``libraries/`` folder, or the model cache in ``~/.cache/atlas-ms``). Copies
in the project folder are deleted with the project. A path typed in the
tables uses the file where it is, without a copy (better for many large raw
files).

Panel's FileDropper keeps each file in memory until it is complete, which a
raw file of several GB would not survive: ``DiskDropper`` writes each piece
to the disk as it arrives instead (a ``.part`` file, renamed when complete).
"""

from pathlib import Path

import panel as pn
import param


class DiskDropper(pn.widgets.FileDropper):
    """A drop area that writes the dropped files into ``folder``; ``saved`` lists them."""

    folder = param.String(default="", doc="Where dropped files are written (created if needed).")
    suffixes = param.List(default=[], item_type=str, doc="Accepted file extensions (lower case); empty: any.")
    saved = param.List(default=[], doc="Paths of the files written, in the order they completed.")
    rejected = param.List(default=[], doc="Names of the dropped files that were refused.")

    # These parameters stay on the server (not sent to the browser's widget).
    _rename = {**pn.widgets.FileDropper._rename, "folder": None, "suffixes": None, "saved": None, "rejected": None}

    def _process_event(self, event) -> None:
        """
        One piece of one file (Panel 1.9 sends an "upload_event" per piece,
        with the file name, the piece number, the number of pieces and the
        bytes). Removing a file from the drop area keeps the copy written.
        """
        if event.event_name != "upload_event":
            return
        data = event.data
        name = Path(data["name"]).name  # a name only: never a folder chosen by the browser
        if not self.folder or (self.suffixes and Path(name).suffix.lower() not in self.suffixes):
            if data["chunk"] == data["total_chunks"]:
                self.rejected = [*self.rejected, name]
            return
        target = Path(self.folder).expanduser() / name
        partial = target.with_name(target.name + ".part")
        target.parent.mkdir(parents=True, exist_ok=True)
        piece = data["data"]
        with open(partial, "wb" if data["chunk"] == 1 else "ab") as handle:
            handle.write(piece.encode() if isinstance(piece, str) else piece)
        if data["chunk"] == data["total_chunks"]:
            partial.replace(target)
            self.saved = [*self.saved, str(target)]


def drop_area(folder: str, suffixes: list[str], label: str, multiple: bool = True) -> DiskDropper:
    """A compact drop area for ``suffixes`` files, written into ``folder``."""
    return DiskDropper(folder=folder, suffixes=suffixes, multiple=multiple, layout="compact",
                       label=label, sizing_mode="stretch_width", height=90)
