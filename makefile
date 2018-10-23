CC = g++
CXXFLAGS = -O3
LDFLAGS = -lm
TARGET = clip
OBJS = clip.o greiner.o polygon.o utilities.o connector.o gpc.o martinez.o

$(TARGET): $(OBJS)

greiner.o: greiner.cpp greiner.h segment.h point.h polygon.h martinez.h utilities.h

polygon.o: polygon.cpp polygon.h segment.h point.h utilities.h

utilities.o: utilities.cpp utilities.h polygon.h segment.h point.h

connector.o: connector.cpp connector.h segment.h point.h martinez.h polygon.h utilities.h

gpc.o: gpc.cpp gpc.h polygon.h segment.h point.h

martinez.o: martinez.cpp martinez.h polygon.h segment.h point.h utilities.h connector.h

clip.o: clip.cpp polygon.h segment.h point.h utilities.h martinez.h connector.h greiner.h gpc.h

clean:
	rm $(TARGET) $(OBJS) *~
