TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

DEFINES += IDE
DEFINES += "FILTER_SIZ=3"

SOURCES += \
        main.c \
    convol.c \
    filter.c \
    matrix.c \
    median.c \
    hist.c \
    pgmfil.c \
    imgoper.c \
    edge.c \
    kmeans.c \
    lst.c \
    arr.c \
    offfil.c \
    proj.c

HEADERS += \
    convol.h \
    filter.h \
    matrix.h \
    utils.h \
    median.h \
    hist.h \
    pgmfil.h \
    imgoper.h \
    edge.h \
    kmeans.h \
    lst.h \
    arr.h \
    offfil.h \
    proj.h
