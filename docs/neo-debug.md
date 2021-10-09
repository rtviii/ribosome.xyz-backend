Know you have a permission problem, your user doesn't have access to some files in the folder /var/lib/neo4j/data/databases/graph.db.

Can you do a sudo chown -R majroud /var/lib/neo4j/data and try again.
Also, be sure that the Neo4j service is not already running.i


https://neo4j.com/developer/kb/permission_denied_errors_after_neo4j_admin/
The recommended approach to all neo4j-admin commands, such as backup, restore, store-info, import is to leverage sudo -u neo4j:

*sudo -u neo4j neo4j-admin backup --from=localhost --name=graph.db_backup_with_user --backup-dir=/tmp*
----

CERTS:https://stackoverflow.com/a/58127788/10697358



