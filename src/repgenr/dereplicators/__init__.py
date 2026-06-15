"""Pluggable dereplication tools (dRep, skDER, galah, sourmash, ...)."""

from .base import Dereplicator, DerepParams, DerepResult, registry

__all__ = ["Dereplicator", "DerepParams", "DerepResult", "registry"]
