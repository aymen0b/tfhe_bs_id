#include <string>
#include "tfhe_bootstrapping.h"

#include <sys/types.h>
#include <unistd.h>

using namespace std;

int main(int argc, char **argv)
{
/*char buf[128];
snprintf(buf, sizeof(buf), "cat /proc/%d/maps", (int)getpid());
system(buf);
*/
  if (argc < 2)
    cout << "[main-error] please enter valid arguments: plain modulus" << endl;
  else {
    int plain_mod = stoi(argv[1], 0, 10);

    cout << "[main-input] plain modulus= " << plain_mod << endl;
    cout << endl;

    //test bootstrapping
    test_bootstrapping_wParams(plain_mod);
    cout << endl;
  }

  return 0;
}
