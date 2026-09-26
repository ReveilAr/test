"""
Annotation: candidate identities for the features, with confidence levels.

Each source writes ``results/annotations/<source>.parquet`` in one shared
format (``schema.py``):

* ``libraries.py``: spectral library search (matchms), your MSP / MGF /
  GNPS JSON libraries;
* ``lipids.py``: rule-based lipid annotation (class fragments and neutral
  losses, species from the precursor m/z, chains when fragments show them);
* (later) MS2Query, SIRIUS.

``harmonize.py`` then combines them: Schymanski confidence levels (it may
lower a source's proposed level, never raise it), flags for disagreements,
the best annotation per feature and a class consensus per molecular family.
"""
