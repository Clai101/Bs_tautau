

#include "sqlite3.h"

int
main() {
  sqlite3_open_v2(":memory:",0,SQLITE_OPEN_READONLY,0);
return 0;
}
