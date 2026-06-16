"""Pluggable SNP typers / variant callers (simple, Snippy, ParSNP, ...)."""

from .base import SnpParams, SnpResult, SnpTyper, registry

__all__ = ["SnpParams", "SnpResult", "SnpTyper", "registry"]
