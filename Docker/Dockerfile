
FROM quay.io/pypa/manylinux2014_x86_64

ENV QT_QPA_PLATFORM minimal

ENV CRYPTOGRAPHY_DONT_BUILD_RUST=1


#_______________________SETTING_PYPI_CREDENTIAL_________________________________
ARG PyPiPassword
ENV PyPiPassword=${PyPiPassword}

ARG PyPiUsername
ENV PyPiUsername=${PyPiUsername}

RUN echo Pip password inputed: ${PyPiPassword}
RUN echo Pip username inputed: $PyPiUsername


#_______________________DOWNLOADING_PACKAGE_____________________________________
RUN mkdir Project
ADD requirements.txt Project/
ADD . Project/SuPyMode
RUN cd Project/SuPyMode && git submodule init && git submodule update
RUN rm -rf Project/SuPyMode/build Project/SuPyMode/dist Project/SuPyMode/extern/eigen/build




#_______________________INSTALL_PYTHON__________________________________________
RUN yum install -y wget cmake nano openssl-devel libffi-devel
RUN cd /usr/src && wget https://www.python.org/ftp/python/3.8.0/Python-3.8.0.tgz && tar xzf Python-3.8.0.tgz
RUN cd /usr/src/Python-3.8.0 && ./configure --enable-optimizations && make altinstall
RUN python3.8 -m pip install --upgrade pip


#_______________________INSTALL_SUPYMODE________________________________________
RUN pip install -r Project/requirements.txt
RUN cd Project/SuPyMode/extern/eigen && mkdir build && cd build && cmake .. && make install
RUN cd Project/SuPyMode/ && mkdir build && cd build && cmake .. && make install



#_______________________UPLOADING_PYPI__________________________________________
RUN echo Compiling and uploading wheel for Python3.8
RUN pip3.8 install -r /Project/SuPyMode/Docker/requirements.txt
RUN cd /Project/SuPyMode/ && python3.8 setup.py bdist_wheel
RUN auditwheel repair /Project/SuPyMode/dist/*.whl -w /Project/SuPyMode/output
RUN make -C /Project/SuPyMode/build UploadPypi
#RUN make -C /Project/SuPyMode CleanProject


#RUN echo Compiling and uploading wheel for Python3.7
#RUN cd /GitProject/PyMieSim && cmake . -DVERSION=3.7 && make clean && make Fibonacci Scatterer _Scatterer _Coupling GaussianBeam
#RUN cd /GitProject/PyMieSim && python3.7 -m pip install -r requirements.txt
#RUN cd /GitProject/PyMieSim/ && python3.7 setup.py bdist_wheel
#RUN auditwheel repair /GitProject/PyMieSim/dist/*37*.whl -w /GitProject/PyMieSim/output
#RUN make -C ./GitProject/PyMieSim UploadPypi
#RUN make -C ./GitProject/PyMieSim CleanProject


#RUN echo Compiling and uploading wheel for Python3.6
#RUN cd /GitProject/PyMieSim && cmake . -DVERSION=3.6 && make clean && make Fibonacci Scatterer _Scatterer _Coupling GaussianBeam
#RUN cd /GitProject/PyMieSim && python3.6 -m pip install -r requirements.txt
#RUN cd /GitProject/PyMieSim/ && python3.6 setup.py bdist_wheel
#RUN auditwheel repair /GitProject/PyMieSim/dist/*36*.whl -w /GitProject/PyMieSim/output
#RUN make -C ./GitProject/PyMieSim UploadPypi
#RUN make -C ./GitProject/PyMieSim CleanProject
