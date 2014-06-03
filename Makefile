JC = javac
JFLAGS = -g
DST = build
SRC = java/lapsolver

all:
	mkdir -p $(DST)
	javac -d $(DST) $(SRC)/*.java
	find $(DST) -name *.class | xargs jar cvf $(DST)/lapsolver.jar

clean:
	$(RM) -rf $(DST)

