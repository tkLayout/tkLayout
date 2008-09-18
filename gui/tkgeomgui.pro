TEMPLATE	= app
LANGUAGE	= C++

CONFIG	+= qt warn_on release

HEADERS	+= gui_constants.h \
	filehandler.h

SOURCES	+= tkgeomgui.cpp \
	FileHandler.cc

FORMS	= maindialog.ui

unix {
  UI_DIR = .ui
  MOC_DIR = .moc
  OBJECTS_DIR = .obj
}

