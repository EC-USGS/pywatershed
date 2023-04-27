class Accessor:
    """A base class for dict access on self."""

    def __setitem__(self, name: str, value) -> None:
        setattr(self, name, value)
        return None

    def __getitem__(self, name: str):
        return getattr(self, name)

    def __delitem__(self, name: str) -> None:
        delattr(self, name)
        return None
