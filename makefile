evolve_A.png : W.py Datos.dat 
	python W.py
    
Datos.dat : W.x
	./W.x

W.x : W.cpp
	c++ W.cpp -o W.x

clean:
	rm W.x *.dat