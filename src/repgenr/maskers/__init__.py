"""Core-SNP alignment maskers (pluggable via the ``repgenr.maskers`` group)."""

from .base import Masker, registry

__all__ = ["Masker", "registry"]
