class DNA:
    __doc__ = """Class representing DNA as a string sequence."""
    basecomplement = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
                      'a':'t', 'c':'g', 'g':'c', 't':'a',
                      'M':'K', 'R':'Y', 'W':'W', 'S':'S', 'Y':'R', 'K':'M',
                      'm':'k', 'r':'y', 'w':'w', 's':'s', 'y':'r', 'k':'m',
                      'V':'B', 'H':'D', 'D':'H', 'B':'V', 
                      'v':'b', 'h':'d', 'd':'h', 'b':'v',
                      'N':'N', 'X':'X',
                      'n':'n', 'x':'x',
                      '-':'-',}
    # IUB encoding:
    # M = A/C, R = A/G, W = A/T, S = C/G, Y = C/T, K = G/T,
    # V = A/C/G, H = A/C/T, D = A/G/T, B = C/G/T, N/X = A/C/G/T, 
    # Celera encoding:
    # m = -/A/C, r = -/A/G, w = -/A/T, s = -/C/G, y = -/C/T, k = -/G/T,
    # v = -/A/C/G, h = -/A/C/T, d = -/A/G/T, b = -/C/G/T, n/x = -/A/C/G/T, 
    
    def __init__(self, s):
        """Create DNA instance initialized to string s."""
        self.seq = s
        return
    def transcribe(self):
        """Return as RNA string."""
        return self.seq.replace('T','U')
    def reverse(self):
        """Return DNA string in reverse order."""
        letters = list(self.seq)
        letters = letters.reverse()
        return ''.join(letters)
    def complement(self):
        """Return the complementary DNA string."""
        letters = list(self.seq)
        letters = [self.basecomplement[base] for base in letters]
        return ''.join(letters)
    def reversecomplement(self):
        """Return the reverse complement of the DNA string."""
        letters = list(self.seq)
        letters.reverse()
        letters = [self.basecomplement[base] for base in letters]
        return ''.join(letters)
    def gc(self):
        """Return the portion of DNA composed of G or C."""
        s = self.seq
        gc = s.count('G') + s.count('C')
        return gc * 1. / len(s)
    def codons(self):
        """Return list of codons for the DNA string."""
        s = self.seq
        end = len(s) - (len(s) % 3) - 1
        codons = [s[i:i+3] for i in range(0,end,3)]
        return codons
