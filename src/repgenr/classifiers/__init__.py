"""Classifier adapters: a genome FASTA in, a GTDB lineage out."""

from .base import Classification, Classifier, ClassifyParams, registry

__all__ = ["Classification", "Classifier", "ClassifyParams", "registry"]
