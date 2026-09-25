"""
Feature-based molecular networking.

    spectra   one MS2 spectrum per feature: loading, cleaning, quality control
    scoring   pairwise similarity (modified cosine or MS2DeepScore) -> candidate edges
    graph     network construction (GNPS algorithm), families, communities, layout, GraphML
"""
