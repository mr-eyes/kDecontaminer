#include "sqliteManager.hpp"
#include "sqlite3pp.h"


int SQLiteManager::callback(void *NotUsed, int argc, char **argv, char **azColName) {
    int i;
    for (i = 0; i < argc; i++) {
        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    }
    printf("\n");
    return 0;
}

void SQLiteManager::create_reads_table() {
    // 1: Single iteration
    // 2: Hierarchical

    string _sqlite_checkTable = "SELECT ID FROM reads LIMIT 1";
    this->rc = sqlite3_exec(this->db.db_, _sqlite_checkTable.c_str(), this->callback, 0, &this->zErrMsg);

    if (this->rc) {

        fprintf(stderr, "`reads` table was not found.\n");
        fprintf(stderr, "Creating `reads` table...\n");

        const char *_sqlite_create_table;
        string _sqlite_create_index;


        _sqlite_create_table = "CREATE TABLE IF NOT EXISTS `reads` ("
                               "`ID`	INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,"
                               "`gene_id`	INTEGER,"
                               "`read_id`	TEXT,"
                               "`PE_seq1`	TEXT,"
                               "`PE_seq2`	TEXT"
                               ");";

        _sqlite_create_index = "CREATE INDEX genes_idx ON reads (gene_id);";

        this->rc = sqlite3_exec(this->db.db_, _sqlite_create_table, this->callback, 0, &this->zErrMsg);
        if (this->rc != SQLITE_OK) {
            fprintf(stderr, "SQL error: %s\n", this->zErrMsg);
            sqlite3_free(this->zErrMsg);
        } else {
            fprintf(stderr, "Table created successfully\n");
        }

        this->rc = sqlite3_exec(this->db.db_, _sqlite_create_index.c_str(), this->callback, 0, &this->zErrMsg);

        if (this->rc != SQLITE_OK) {
            fprintf(stderr, "SQL error: %s\n", this->zErrMsg);
            sqlite3_free(this->zErrMsg);
        } else {
            fprintf(stderr, "DB Index created successfully.\n");
        }

    } else {
        fprintf(stderr, "`reads` table found.\n");
    }

    fprintf(stderr, "Done initializing DB.\n");
    string _sqlite_sync = "PRAGMA synchronous = OFF;";
    this->rc = sqlite3_exec(this->db.db_, _sqlite_sync.c_str(), this->callback, 0, &this->zErrMsg);
    if (this->rc != SQLITE_OK) {
        fprintf(stderr, "SQL error: %s\n", this->zErrMsg);
        sqlite3_free(this->zErrMsg);
    }

}

bool SQLiteManager::check_reads_table() {
    string _sqlite_checkTable = "SELECT ID FROM reads LIMIT 1";
    this->rc = sqlite3_exec(this->db.db_, _sqlite_checkTable.c_str(), this->callback, nullptr, &this->zErrMsg);
    return this->rc == 0;
}

SQLiteManager::SQLiteManager(const string &db_file) {
    this->db = sqlite3pp::database(db_file.c_str());
}

void SQLiteManager::close() {
    sqlite3_close(this->db.db_);
}