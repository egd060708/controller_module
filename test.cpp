#include "mpcCal.h"
#include "PIDmethod.h"

int main()
{
    mpcCal<6, 12, 20, 1, 5> balanceController;
    PIDmethod pid;
    return 0;
}