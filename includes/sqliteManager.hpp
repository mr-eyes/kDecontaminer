#include <stdio.h>
#include <sqlite3.h>
#include <string>
#include "iostream"
#include "sqlite3pp.h"

using namespace std;

class SQLiteManager {

public:
    char *zErrMsg = 0;
    int rc;

    static int callback(void *NotUsed, int argc, char **argv, char **azColName);


public:
    sqlite3pp::database db;
    SQLiteManager(const string& db_file);
    void create_reads_table();
    bool check_reads_table();
    void close();

};