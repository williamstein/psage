"""
This module implements a compressed key to Python object store using MongoDB.

This module defines one class ObjectDB, which we instantiate using a
MongoDB database db.  The resulting instance then works somewhat like
a dictionary, except that objects are pickled and stored on disk in
MongoDB.
"""

class ObjectDB:
    def __init__(self, db):
        from gridfs import GridFS
        self.gridfs = GridFS(db)

    def __setitem__(self, key, obj):
        self.save(obj, key)

    def __getitem__(self, key):
        return self.load(key)

    def __delitem__(self, key):
        from pymongo.objectid import ObjectId
        if not isinstance(key, ObjectId):
            id = self.gridfs.get_last_version(key)._id
        else:
            id = key
        self.gridfs.delete(id)

    def __repr__(self):
        return "Key-value database"

    def keys(self):
        """Return list of filenames of objects in the gridfs store."""
        return self.gridfs.list()

    def object_ids(self):
        """Return list of id's of objects in the gridfs store, which
        are not id's of objects with filenames."""
        v = self.gridfs._GridFS__files.find({'filename':{'$exists':False}},['_id'])
        return [x['_id'] for x in v]

    def has_key(self, key):
        return self.gridfs.exists(filename=key)

    def save(self, obj, key=None, compress=None):
        """Save Python object obj to the grid file system self.gridfs.
        If key is None, the file is stored by MongoDB assigned
        ObjectID, and that id is returned.
        """
        from sage.all import dumps
        data = dumps(obj, compress=compress)
        if key is not None:
            self.gridfs.put(data, filename=key)
            return key
        else:
            # store by MongoDB assigned _id only, and return that id.
            return self.gridfs.put(data)

    def load(self, key, compress=True):
        from pymongo.objectid import ObjectId
        if isinstance(key, ObjectId):
            data = self.gridfs.get(key).read()
        else:
            data = self.gridfs.get_last_version(key).read()
        from sage.all import loads
        return loads(data, compress=compress)
