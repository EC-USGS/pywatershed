class Accessor:
    """A base class for dict access on self.

    This class is used by most other classes to ``__setitem__`,
    ``__get_item__``, and ``__delitem__`` in a ``dict`` style

    """

    def __setitem__(self, name: str, value) -> None:
        setattr(self, name, value)
        return None

    def __getitem__(self, name: str):
        return getattr(self, name)

    def __delitem__(self, name: str) -> None:
        delattr(self, name)
        return None
