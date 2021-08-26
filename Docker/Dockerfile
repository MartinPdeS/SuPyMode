
FROM quay.io/pypa/manylinux2014_x86_64

ENV QT_QPA_PLATFORM minimal

ENV CRYPTOGRAPHY_DONT_BUILD_RUST=1



ARG PyPiPassword
ENV PyPiPassword=${PyPiPassword}

ARG PyPiUsername
ENV PyPiUsername=${PyPiUsername}

RUN echo Pip password inputed: ${PyPiPassword}
RUN echo Pip username inputed: $PyPiUsername

RUN mkdir Project

ADD requirements.txt Project/
RUN cd Project && git clone https://github.com/MartinPdeS/SuPyMode.git



#RUN yum install -y wget cmake nano openssl-devel libffi-devel

#RUN cd /usr/src && wget https://www.python.org/ftp/python/3.8.0/Python-3.8.0.tgz && tar xzf Python-3.8.0.tgz

#RUN cd /usr/src/Python-3.8.0 && ./configure --enable-optimizations && make altinstall

#RUN python3.8 -m pip install --upgrade pip

#RUN pip install -r requirements.txt


#RUN echo Copying and cleaning files...
#RUN mkdir GitProject && mkdir GitProject/SuPyMode
#ADD . ./GitProject/SuPyMode/

#RUN cd GitProject/SuPyMode/ && cmake . && cd ../..
#RUN make -C ./GitProject/SuPyMode CleanProject


#RUN echo Installing Boost dependencies...
#RUN cp -r /GitProject/PyMieSim/extern/math/include/boost /usr/include


#RUN echo Installing Fortran 90 dependencies...
#RUN cmake ./GitProject/PyMieSim && make -C ./GitProject/PyMieSim installBessel


#RUN echo Compiling and uploading wheel for Python3.8
#RUN cd /GitProject/PyMieSim && cmake . -DVERSION=3.8 && make clean && make Fibonacci Scatterer _Scatterer _Coupling GaussianBeam
#RUN cd /GitProject/PyMieSim && python3.8 -m pip install -r requirements.txt
#RUN cd /GitProject/PyMieSim/ && python3.8 setup.py bdist_wheel
#RUN auditwheel repair /GitProject/PyMieSim/dist/*38*.whl -w /GitProject/PyMieSim/output
#RUN make -C ./GitProject/PyMieSim UploadPypi
#RUN make -C ./GitProject/PyMieSim CleanProject


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
