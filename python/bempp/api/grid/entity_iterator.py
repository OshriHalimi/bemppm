"""Define the EntityIterator class."""

from .entity import Entity as _Entity


class EntityIterator(object):
    """Implements iterators over entities."""

    def __init__(self, codimension, impl):
        self._codimension = codimension
        self._impl = impl

    def next(self):
        return _Entity(self._codimension,
                       self._impl.__next__())

    def __next__(self):
        return self.next()

    def __iter__(self):
        return self
