#include <qapplication.h>
#include "maindialog.h"

int main( int argc, char ** argv )
{
    QApplication a( argc, argv );
    MainDialog w;
    w.show();
    a.connect( &a, SIGNAL( lastWindowClosed() ), &a, SLOT( quit() ) );
    return a.exec();
}
