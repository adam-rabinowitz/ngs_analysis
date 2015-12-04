""" This module is built to handle the output of functions. The use of
the 'handle_output' function within a function allows that function to
return output as a list or to a file or down a multiprocessing pipe.
"""


class OutputList(object):
    """ Generates an object to handle the output of function as a list.
    Object is initialised with a list. The two defined functions 'add'
    and 'close' add elements to the list and return the list
    respectively.
    """
    
    def __init__(self):
        self.output = list()
    
    def add(self, data):
        self.output.append(data)
    
    def close(self):
        return self.output


class OutputFile(object):
    """ Generates an object to handle the output of function to a file.
    Object is initialised by opening a file using the supplied file
    name. If the file name ends with '.gz' then it is assumed the file
    is gzipped. The function 'add' writes a supplied object as a line in
    the file. If the object is a list or tuple then the the indvidual
    elements are converted to strings and joined with an intervening
    tab. Otherwise, the object is converted to a string and added to
    the file. The function 'close' closes the file and returns the file
    handle.
    """
    
    def __init__(self, output):
        if output.endswith('.gz'):
            import gzip
            self.output = gzip.open(output, 'w')
        else:
            self.output = open(output, 'w')
    
    def add(self, data):
        if isinstance(data, list) or isinstance(data, tuple):
            data = map(str, data)
            self.output.write('\t'.join(data) + '\n')
        else:
            data = str(data)
            self.output.write(data + '\n')
    
    def close(self):
        self.output.close()
        return self.output


class OutputPipe(object):
    """ Generates a class to handle output of a function to a pipe.
    Object is initialised with an open pipe. The function 'add' sends
    an object down the pipe. The function 'close' closes the pipe and
    returns it.
    """
    
    def __init__(self, output):
        self.output = output
    
    def add(self, data):
        self.output.send(data)
    
    def close(self):
        self.output.close()
        return self.output


class InputList(object):
    """ Generates an object to handle using a list as input to a
    function. Object is initialised with a list. The two defined
    'next' and 'close' return elements from the list and return the
    list respectively. When the list is empty an EOFerror is raised.
    """
    
    def __init__(self, listin):
        self.input = iter(listin)
    
    def next(self):
        try:
            return next(self.input)
        except StopIteration:
            raise EOFError('End Of List')
    
    def close(self):
        return self.input


class InputFile(object):
    """ Generates an object to handle the input of a file to a
    function. Object is initialised by opening a file using a supplied
    file name. If the file name ends with '.gz' it is assumed the file
    is gzipped. The function 'next' sequentially extracts a line from
    the file and removes newline characters. Line containing tabs are
    split into a list using the tabs as delimiters. If the line is
    blank the file is closed and an EOFError is raised. The function
    'close' closes the file and returns the file handle.
    """
    
    def __init__(self, filein):
        if filein.endswith('.gz'):
            import gzip
            self.input = gzip.open(filein, 'r')
        else:
            self.input = open(filein, 'r')
    
    def next(self):
        data = self.input.readline().strip()
        if data:
            data = data.split('\t')
            if len(data) == 1:
                return data[0]
            else:
                return data
        else:
            self.input.close()
            raise EOFError('End Of File')
    
    def close(self):
        self.input.close()
        return self.input


class InputPipe(object):
    """ Generates a class to allow a pipe to be used as input to a
    function. Object is initialised with an open pipe. The function
    'next' extracts the next object from the pipe. If the pipe is empty
    and the other end is open then the function will wait until the an
    object is sent. If the pipe is empty and the other end is closed
    then and EOF error is raised. The function 'close' closes the pipe
    and returns it.
    """
    
    def __init__(self, pipein):
        self.input = pipein
    
    def next(self):
        data = self.input.recv()
        if data == None:
            raise EOFError('No More Data In Pipe')
        else:
            return data
    
    def close(self):
        self.input.close()
        return self.input


def handleout(Output):
    """ Function takes a single argument termed 'output'. If 'output'
    is a string then a file output object (OutputFile) is returned. If
    the 'output' is a list then a list output object (OutputList) is
    returned. If 'output' is a pipe then a pipe output object
    (OutputPipe) is returned. All objects have the functions 'add' and
    'close'.The function 'add' adds another item to the output while
    'close' closes the file/pipe and returns the file handle, pipe, or
    list.
    """
    if isinstance(Output, str):
        outObject = OutputFile(Output)
    elif isinstance(Output, list):
        outObject = OutputList()
    else:
        from _multiprocessing import Connection
        if isinstance(Output, Connection):
            outObject = OutputPipe(Output)
        else:
            raise TypeError('Output must be a file-name, list or pipe')
    return outObject


def handlein(Input):
    """ Function takes a single argument termed 'input'. If 'input'
    is a string then a file input object (InputFile) is returned. If
    the 'input' is a list then a list input object (InputList) is
    returned. If 'input' is a pipe then a pipe input object
    (InputPipe) is returned. All objects have the functions 'next' and
    'close'. The 'next' function returns the next input item and raises
    the EOFError exception when there is no more data available. The
    function 'close' closes the file/pipe and returns the file handle, pipe
    or list.
    """
    if isinstance(Input, str):
        inObject = InputFile(Input)
    elif isinstance(Input, list):
        inObject = InputList(Input)
    else:
        from _multiprocessing import Connection
        if isinstance(Input, Connection):
            inObject = InputPipe(Input)
        else:
            raise TypeError('Input must be a file-name, list or pipe')
    return inObject
