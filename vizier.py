import os
import numpy as np

class catalog():
    def __init__(self, directory, onlyFiles=None, verbose=False):
        self.directory = directory
        # -- readme:
        try:
            self._readme = open(os.path.join(directory, 'ReadMe')).readlines()
        except:
            raise NameError(os.path.join(directory, 'ReadMe')+' does not exist')
        self.rawdata = {}
        self.data = {}
        self.verbose = verbose
        self._readReadMe()
        if onlyFiles is None:
            # -- default is all files
            onlyFiles = list(self.decode.keys())
        for f in onlyFiles:
            try:
                print('reading', f)
                self._loadFile(f)
                self._reformatData(f)
            except:
                pass

    def _readReadMe(self):
        """
        sets the decoding dictionnary
        """
        self.decode = {}
        fileName = None
        k = 0
        if self.verbose:
            print(len(self._readme), 'lines in readme')
        while k < len(self._readme):
            #if self.verbose:
            #    print k, self._readme[k]
            if 'Description of file:' in self._readme[k]:
                fileName = self._readme[k].split(':')[1].strip()
                if self.verbose:
                    print('         adding file', fileName)
                k += 4
                tmp = []
            if not fileName is None:
                l = self._readme[k]
                if l[0] == '-':
                    self.decode[fileName] = tmp
                    fileName = None
                elif l[:15].strip()!='':
                    # i_min, i_max, type, unit
                    if not '-' in l[:10]:
                        imin = int(l.split()[0])-1
                        imax = imin+1
                    else:
                        imin = int(l[:4])-1
                        imax = int(l[5:8])
                    tmp.append(( imin, imax,
                                l[10:15].strip(), l[16:22].strip(),# type, unit
                                l[22:].split()[0].strip(), # name
                                l[33:-1])) # comment
            k += 1
        return

    def _loadFile(self, filename):
        """
        basic loading
        """
        if filename not in self.decode:
            raise ErrorName('unknown file: '+filename)
        lines = open(os.path.join(self.directory, filename)).readlines()
        data = []
        verbose=False
        for l in lines:
            data.append(stripLine(l, self.decode[filename],
                                  verbose=verbose))
            verbose = False
        self.rawdata[filename] = data
        return

    def _reformatData(self, filename):
        """
        """
        if filename not in self.decode:
            raise ErrorName('unknown file: '+filename)
        print('reformating', filename)
        data = {}
        for k,f in enumerate(self.decode[filename]):
            tmp = [x[k] for x in self.rawdata[filename]]
            if 'I' in f[2] or 'F' in f[2]:
                data[f[4]] = np.array(tmp)
            else:
                data[f[4]] = tmp
        # check if coordinates are present, if so create RA.h and DE.d
        if all([k in data for k in ['RAh', 'RAm', 'RAs',
                                          'DE-','DEd', 'DEm', 'DEs']]):
            data['RA.h'] = data['RAh']+data['RAm']/60.0+data['RAs']/3600.0
            data['DEC.d'] = data['DEd']+data['DEm']/60.0+data['DEs']/3600.0
            data['DEC.d'][np.array(data['DE-'])=='-'] *=-1

        self.data[filename] = data
        self.rawdata.pop(filename)
        return



def stripLine(l, form, verbose=False):
    tmp = []
    # slice
    if verbose:
        print(l)
    for f in form:
        tmp.append(l[f[0]:f[1]])
        if verbose:
            print(' '*3, f, tmp[-1], '(', len(tmp[-1]), ')')
        if 'I' in f[2]:
            try:
                tmp[-1] = int(tmp[-1])
            except:
                tmp[-1] = np.nan

        if 'F' in f[2]:
            try:
                tmp[-1] = float(tmp[-1])
            except:
                tmp[-1] = np.nan
    if verbose:
        print('')
        print(tmp)
    return tmp