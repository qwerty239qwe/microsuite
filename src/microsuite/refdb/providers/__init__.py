from __future__ import annotations

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers._base import RefDbProvider

_PROVIDERS: dict[str, RefDbProvider] = {}


def register_provider(provider: RefDbProvider) -> None:
    _PROVIDERS[provider.name] = provider


def get_provider(name: str) -> RefDbProvider:
    provider = _PROVIDERS.get(name)
    if provider is None:
        choices = ", ".join(sorted(_PROVIDERS)) or "(none registered)"
        raise MicrobiomeSuiteError(f"Unknown reference-DB provider '{name}'. Available: {choices}")
    return provider
