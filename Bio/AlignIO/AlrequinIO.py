# Copyright

"""
AlignIO support for the "Arlequin" format used in Laurent Excoffier's population
genetic package of the same name.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).
"""
#import string

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Interfaces import AlignmentIterator, SequentialAlignmentWriter

try:
    any
except NameError:
    #Hack for Python 2.4
    def any(iterable):
        for element in iterable:
            if element:
               return True
        return False




class ArlequinWriter(SequentialAlignmentWriter):
    """Arlequin alignment writer."""

    def write_alignment(self, alignment):
        """
        """
        pass
        
class ArlequinIterator(AlignmentIterator):
    """Reads a Arlequin alignment file returning a MultipleSeqAlignment iterator.

    At present 
    """

    def _clean(self, line):
        """Remove comments from a line """
        return line.split('#')[-1].strip()

    def _get_value(self, line):
        """get attributes from the header sections of file
        """
        return line.split('=')[1].replace('"', "").replace("'", "")

    def _is_profile(self):
        line = self.handle.readline()
        while '[Profile]' not in self._clean(line):
            if not line:
                return 
            line = self.handle.readline()
        return True

    def _get_profile(self):
        line = self.handle.readline()
        while '[Data]' not in self._clean(line):
            line = self._clean(self.handle.readline())
            if 'DataType' in line:
                dtype = self._get_value(line)
                if dtype != 'DNA':
                    raise ValueError(
                    'At the moment arlequin parser only supports haplotypic datafiles this file contains {0} data'.format(dtype))
            if 'MissingData' in line:
                self.alphabet = Gapped(IUPAC.ambiguous_dna, self._get_value(line))
            if "NbSamples" in line:
                try:
                    self.nsamples = int(self._get_value(line))
                except ValueError:
                   print 'Cannot coerce NbSamples from profile to an integer'
         #TODO
         #needs something to set alphabet if missing data is missing
        
    def next(self):
        handle = self.handle
        line = handle.readline()
        self._is_profile()
        #we have hit a profile, so lets read it
        self._get_profile()
        #
        alignment = MultipleSeqAlignment([], self.alphabet)
        for sample in range(0, self.nsamples):
            line = self._clean(self.handle.readline())
            while "SampleName" not in line:
                line = self._clean(self.handle.readline())
            sample_name = self._get_value(line)
            line = self._clean(self.handle.readline())
            n_ind = int(self._get_value(line))
            line = self._clean(self.handle.readline())
            n = 0
            while n < n_ind:
                line = self._clean(self.handle.readline())
                if line: #skip empty lines
                    id, freq, seq = line.split()
                    d =  "freq={0} sample={1}".format(freq, sample_name)
                    alignment.append(SeqRecord(seq=Seq(seq, self.alphabet), \
                                               id = id, name = id,          \
                                               description = d))
                    record = alignment[-1]
                    record.annotations["sample"] = sample_name
                    record.annotations["frequency"] = int(freq)
                    n += int(freq)
        return alignment
        
