from typing import Union

class BIGSDbDatabaseAPIException(Exception):
    pass


class NoBIGSdbMatchesException(BIGSDbDatabaseAPIException):
    def __init__(self, database_name: str, database_schema_id: int, *args):
        super().__init__(f"No matches found with schema with ID {database_schema_id}  in the database \"{database_name}\".", *args)

class NoBIGSdbExactMatchesException(NoBIGSdbMatchesException):
    def __init__(self, database_name: str, database_schema_id: int, *args):
        super().__init__(f"No exact match found with schema with ID {database_schema_id}  in the database \"{database_name}\".", *args)

class NoSuchBIGSdbDatabaseException(BIGSDbDatabaseAPIException):
    def __init__(self, database_name: str, *args):
        super().__init__(f"No database \"{database_name}\" found.", *args)

class NoSuchBigSdbSchemaException(BIGSDbDatabaseAPIException):
    def __init__(self, database_name: str, database_schema_id: int, *args):
        super().__init__(f"No schema with ID {database_schema_id}  in \"{database_name}\" found.", *args)
