update-cloud-db.sh lets you specify a database and an optional list of collections that you with to copy from m0 to ms.

By default, it will look in db-collections.txt for the specified database and update all collections listed there. If you manually specify a collection it will verify that the collection is listed in db-collections.txt for the specified database (this is primarily to protect us from making typos).
