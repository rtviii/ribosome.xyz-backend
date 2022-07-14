


```zsh
neoconf ="sudo vim /etc/neo4j/neo4j.conf"
neolog  ="sudo vim /var/log/neo4j"
neodata ="cd /var/lib/neo4j/data"
neologin="sudo -u neo4j /bin/cypher-shell -u neo4j -p 55288"
```


CERTS: https://stackoverflow.com/a/58127788/10697358


- get apoc : [ `https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases` ](https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases) <--from here. `wget` the "all" `.jar` into `/var/lib/neo4j/plugins`. Restart the database.

- `neo4j.conf` example is provided:
  - make sure the default db is always up
  - memory heap is `1024m` or upwards
  - `dbms.default_list_address` is `0.0.0.0`

  Enable file import:

  ```conf
	apoc.import.file.enabled=true
	apoc.import.file.use_neo4j_config=true
  ```

  Apoc and other procedures:

  ```conf
	# full access to the database through unsupported/insecure internal APIs.
	dbms.security.procedures.unrestricted=apoc.*
	# A comma separated list of procedures to be loaded by default.
	# Leaving this unconfigured will load all procedures found.
	dbms.security.procedures.allowlist=apoc.coll.*,apoc.load.*,gds.*,apoc.*
  ```

- change password upon initial login
- create a new database
- move resources into /var/lib/neo4j/import
- create ontologies, constraints
- upload structures










Know you have a permission problem, your user doesn't have access to some files in the folder /var/lib/neo4j/data/databases/graph.db.

Can you do a sudo chown -R majroud /var/lib/neo4j/data and try again.
Also, be sure that the Neo4j service is not already running.i


https://neo4j.com/developer/kb/permission_denied_errors_after_neo4j_admin/
The recommended approach to all neo4j-admin commands, such as backup, restore, store-info, import is to leverage sudo -u neo4j:

*sudo -u neo4j neo4j-admin backup --from=localhost --name=graph.db_backup_with_user --backup-dir=/tmp*
----

CERTS:https://stackoverflow.com/a/58127788/10697358


