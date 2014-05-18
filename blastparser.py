#!/usr/bin/env python
###############################################################################
#
# __script_name__.py - description!
#
###############################################################################
# #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or #
# (at your option) any later version. #
# #
# This program is distributed in the hope that it will be useful, #
# but WITHOUT ANY WARRANTY; without even the implied warranty of #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the #
# GNU General Public License for more details. #
# #
# You should have received a copy of the GNU General Public License #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
# #
###############################################################################

__author__ = "Josh Daly"
__copyright__ = "Copyright 2014"
__credits__ = ["Josh Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = ""
__status__ = "Development"

###############################################################################

import argparse
import sys

from multiprocessing import Pool
from subprocess import Popen, PIPE

from pyparsing import *
from cStringIO import StringIO
import re

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

  # classes here
class BlastHitSummary:
    """
    BlastHitSummary.

    A class used to store information from the BLAST summary header
    (basically, abbreviated sequence name, total score, and total
    expectation).
    """
    __slots__ = ['subject_name', 'total_score', 'total_expect', 'N']
    def __init__(self, subject_name, total_score, total_expect, N):
        self.subject_name = subject_name
        self.total_score = int(total_score)
        self.total_expect = float(total_expect)
        self.N = N

    def __repr__(self):
        return "<BlastHitSummary(%s, total_expect=%g, total_score=%d)>" % \
               (self.subject_name, self.total_expect, self.total_score)
    

class BlastSubjectSubmatch:
    """
    BlastSubjectSubmatch.
    
    A specific submatch (score/alignment) of a query sequence to a
    subject sequence.

    Attributes:

     - expect
     - frame1
     - frame2
     - bits
     - bits_max
     - query_start
     - query_end
     - subject_start
     - subject_end
     - query_sequence
     - subject_sequence
     - alignment

    Usage: ::

        print submatch_obj.expect
        
    (etc.)
     
    """
    __slots__ = ['expect', 'frame1', 'frame2', 'bits', 'bits_max',
                 'query_start', 'query_end', 'query_sequence'
                 'subject_start', 'subject_end', 'subject_sequence',
                 'alignment']
    
    def __init__(self, expect, frame1, frame2, bits, bits_max,
                 q_start, q_end, q_seq, s_start, s_end, s_seq, align):
        self.expect = float(expect)
        self.frame1 = frame1
        self.frame2 = frame2
        self.bits = int(bits)
        self.bits_max = int(bits_max)
        self.query_start = q_start
        self.query_end = q_end
        self.query_sequence = q_seq

        self.subject_start = s_start
        self.subject_end = s_end
        self.subject_sequence = s_seq

        self.alignment = align

    def __repr__(self):
        return "<BlastSubjectSubmatch(expect=%g, query %d-%d, subject %d-%d))>"\
               % (self.expect, self.query_start, self.query_end,
                self.subject_start, self.subject_end)

class BlastSubjectHits:
    """
    BlastSubjectHits.

    A list of all of the matches between a query sequence and a subject
    sequence.

    Attributes:
     * subject_name -- name of subject sequence.
     * subject_length -- length of subject sequence.
     * total_expect -- total expected value of this match
     * total_score -- total score
     * matches -- list of BlastSubjectSubmatch objects.

    Usage: ::

        print hits_object.subject_name, hits_object.subject_length
        for match in hits_object:
           print match
    """
    __slots__ = ['subject_name', 'subject_length', 'matches', 'total_expect',
                 'total_score']
    def __init__(self, subject_name, subject_length, matches):
        self.subject_name = str(subject_name)
        self.subject_length = int(subject_length)
        self.matches = matches

    def __getitem__(self, i):
        return self.matches[i]

    def __len__(self):
        return len(self.matches)

    def __repr__(self):
        seqname = build_short_sequence_name(self.subject_name)
        return "<BlastSubjectHits(%s, %s, %d matches)>" % (seqname,
                                                           self.subject_length,
                                                           len(self))

class BlastDatabase:
    """
    Information about the BLAST database used in a particular Query.
    
    Attributes:
    
      * name -- name of the BLAST database, as given by BLAST.
      * num_seqs -- number of sequences in the database.
      * total_chars -- total number of characters in all sequences combined.
      
    """
    __slots__ = [ 'name', 'num_seqs', 'total_chars' ]
    def __init__(self, name, num_seqs, total_chars):
        """
        BlastDatabase(name, num_seqs, total_chars)
        """
        if not name:
            raise Exception("no database name given")
        self.name = name
        self.num_seqs = int(num_seqs)
        self.total_chars = int(total_chars)

    def __repr__(self):
        return "< BlastDatabase(%s) >" % (self.name,)

class BlastQuery:
    """
    A BLAST query (single sequence against database) containing all results.
    
    Attributes:

      * program -- the BLAST program used (BLASTX, etc.)
      * version -- the version of the BLASTALL binary
      * date -- the date of the BLASTALL binary
      * query_name -- name of query sequence (following 'Query=').
      * query_length -- total length of query sequence.
      * database -- a BlastDatabase instance; has name, # seqs, # characters.
      * hits -- a list of BlastSubjectHits, containing each match + alignment.
      
    Usage: ::

        print query_object.query_name, query_object.query_length
        for hits_object in query_object:
           print hits_object.subject_name
    """
    __slots__ = ['query_name', 'query_length', 'database', '_hitlines','hits',
                 'program', 'date', 'version']
    def __init__(self, program, date, version,
                 query_name, query_length, database, hitlines, hits):
        self.program = program
        if not query_name:
            raise Exception("no query sequence name given")
        self.query_name = query_name
        self.query_length = int(query_length)
        if not isinstance(database, BlastDatabase):
            if database:
                raise Exception("incorrect database object type")
            else:
                raise Exception("no BLAST database information available")
        self.database = database

        self._hitlines = list(hitlines)
        self.hits = list(hits)

    def __repr__(self):
        query_short = build_short_sequence_name(self.query_name)
        return "<BlastQuery(%s, %s, %s)>" % (query_short,
                                             self.query_length,
                                             self.database.name)

    def __len__(self):
        return len(self.hits)

    def __getitem__(self, i):
        return self.hits[i]

class BlastFile:
    CHUNKSIZE=10*1024*1024

    blast_marker = re.compile('^(BLASTX|BLASTN|BLASTP|TBLASTX|TBLASTN)',
                              re.MULTILINE)
    query_marker = '\nQuery='
    no_hits_marker = "***** No hits found ******"
    
    def __init__(self, fp, ignore_no_hits=False):
        self.first = True
        self.fp = fp
        self.ignore_no_hits = ignore_no_hits
        self.total = 0

        self.s = ""
        self.pos = 0
        self.stopped = False

        self._more()

        if not self.blast_marker.match(self.s):
            raise Exception("this doesn't look like a BLAST file")

    def _more(self):
        r = self.fp.read(self.CHUNKSIZE)
        if not r:
            raise StopIteration
        self.s = self.s[self.pos:] + r
        self.pos = 0

    def __iter__(self):
        return self

    def next(self):
        if self.stopped:
            raise StopIteration
            
        while 1:
            a = self.s.find(self.query_marker, self.pos)
            if a >= 0:
                b = self.s.find(self.query_marker, a + 1)

            # do we have a complete record?
            if a >= 0 and b > 0:
                record = self.s[self.pos:b]
                self.pos = b
                self.total += 1

                if self.ignore_no_hits and \
                       record.find(self.no_hits_marker) >= 0:
                    if self.first:
                        self.first = False
                        return record # get database info
                else:
                    return record
                
            # no, we need to get some more from the file.
            else:
                try:
                    self._more()
                except StopIteration:
                    self.stopped = True
                    record, self.s = self.s[self.pos:], ""
                    self.pos = 0
                    return record

class _BlastShelf:
    def __init__(self, filename, mode='r'):
        from shelve import BsdDbShelf
        from bsddb import btopen

        _db = btopen(filename, 'r')
        self.db = BsdDbShelf(_db)

    def __iter__(self):
        db = self.db
        last_k, _ = db.last()
        k, v = db.first()
        while k != last_k:
            yield k, v
            k, v = db.next()
        yield k, v

def open_shelf(filename, mode='r'):
    from shelve import BsdDbShelf
    from bsddb import btopen

    return _BlastShelf(filename, mode)

def parse_file(filename, **kw):
    """
    Parse records from a given file; 'filename' is the path to the file.
    """
    b = BlastParser()
    for record in b.parse_file(filename, **kw):
        yield record

def parse_fp(fp, **kw):
    """
    Parse records out of the given file handle.
    """
    b = BlastParser()
    
    for record in b.parse_fp(fp, **kw):
        yield record

def parse_string(s, **kw):
    """
    Parse records out of a string buffer.
    """
    b = BlastParser()

    for record in b.parse_string(s, **kw):
        yield record
    
class BlastParser:
    """
    BlastParser objects coordinate the use of pyparsing parsers to
    parse complete BLAST records.

    Attributes:

      * blast_record -- an individual BLAST record; returns BlastQuery object.
      * blast_output -- list of records; returns list of BlastQuery objects.

    Methods:

      * reset() -- clear the blast parser of persistent information.
      * parse_string(s)
      * parse_file(filename)
      * parse_fp(fp)
     """
    def __init__(self):
        self.reset()
        self._construct_parser()

    def reset(self):
        self.blast_database = None

    def parse_string(self, s, **kw):
        fp = StringIO(s)
        return self.parse_fp(fp, **kw)

    def parse_file(self, filename, **kw):
        fp = open(filename)
        return self.parse_fp(fp, **kw)

    def parse_fp(self, fp, ignore_errors=True, ignore_no_hits=False):
        self.reset()

        for record in BlastFile(fp, ignore_no_hits):
            errfp = None
            try:
                parse_results = self.blast_record.parseString(record)
            except ParseException:
                print 'ERROR', record[0:50]
                if ignore_errors:

                    if errfp == None:
                        errfp = open('err.txt', 'w')
                    print >>errfp, record

                    continue
                
                raise
                
            yield parse_results[0]

        if errfp:
            errfp.close()

    def set_blast_database(self, x):
        dbname = x.database_name[0].strip()

        database = BlastDatabase(dbname, x.database_num_seqs,
                                 x.database_total_seq)
        self.blast_database = database
        assert self.blast_database

    def set_blast_info(self, x):
        self.blast_program = x.blast_type
        self.blast_version = x.blast_version
        self.blast_date = x.blast_date

    def _construct_parser(self):
        #### BLAST HEADER ####

        # parse:
        """
        BLASTX 2.2.14 [May-07-2006]
        """

        self.program_line = oneOf("BLASTX BLASTP BLASTN TBLASTX TBLASTN" \
                                  ).setResultsName("blast_type") + \
                        Word(nums + ".").setResultsName("blast_version") + \
                        Word("[-]" + alphanums).setResultsName("blast_date")
        self.program_line.setParseAction(self.set_blast_info)

        blast_reference = (Literal("""\
Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, 
Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), 
"Gapped BLAST and PSI-BLAST: a new generation of protein database search
programs",  Nucleic Acids Res. 25:3389-3402.
""") + Optional("""\
Reference for compositional score matrix adjustment: Altschul, Stephen F., 
John C. Wootton, E. Michael Gertz, Richa Agarwala, Aleksandr Morgulis,
Alejandro A. Schaffer, and Yi-Kuo Yu (2005) "Protein database searches
using compositionally adjusted substitution matrices", FEBS J. 272:5101-5109.
""")).suppress()

        # parse:
        """
        Database: nr-divvy.Archaea.seq 
        107,277 sequences; 30,770,611 total letters
        """

        self.db_info_line = \
                     db_line = Literal("Database: ") + \
                     SkipTo(LineEnd() + named_comma_int('') + "sequences;",
                            include=False).setResultsName('database_name') + \
                            named_comma_int('database_num_seqs') + \
                            "sequences;" + \
                            named_comma_int('database_total_seq') + \
                            "total letters"

        # parse:
        """
        Searching..................................................done
        """

        self.searching = (Literal("Searching") + ZeroOrMore('.')).suppress() +\
                         ZeroOrMore(Literal("done"))

        self.blast_top = self.program_line + \
                         blast_reference

        self.db_info = self.db_info_line + self.searching
        self.db_info.setParseAction(self.set_blast_database)

        # parse:
        """
        Query= EKJFPUJ01A000T length=117 xy=0304_2523 region=1
        run=R_2007_02_22_16_50_08_
        (117 letters)
        """

        query_size = LineEnd() + "(" + named_comma_int('query_size') + \
                     "letters)"
        query_name = SkipTo(query_size).setResultsName('query_name')

        query_line = Literal("Query=") + query_name + query_size 
        self.query_line = query_line.setResultsName('query_header')
        

        # parse:
        """
        ***** No hits found ******
        """

        self.no_hits = Literal("***** No hits found ******").suppress()
        self.no_hits = self.no_hits.setResultsName('hits')

        # parse:
        """
                                                                 Score    E
Sequences producing significant alignments:                      (bits) Value
        
ref|YP_136807.1| hypothetical protein rrnAC2273 [Haloarcula mari...    36   0.011
ref|YP_134259.1| hypothetical protein pNG5141 [Haloarcula marism...    28   3.0  
ref|YP_134121.1| hypothetical protein pNG4041 [Haloarcula marism...    27   5.1  
gb|ABI23033.1| succinyl-diaminopimelate desuccinylase [unculture...    27   8.7  
        """

        hits_header = SkipTo((Literal("Value") + NotAny(LineEnd()) + "N") | "Value", include=True)\
                      + OneOrMore(LineEnd())

        endlinepattern = NotAny(LineEnd()) + \
     (named_int('total_score') + named_float('total_expect') + named_int('N') |
      named_int('total_score') + named_float('total_expect')) + LineEnd()

        hitline = NotAny(LineEnd()) + SkipTo(endlinepattern, include=True)
        hitline = hitline.setParseAction(process_hitline)
        
        hitlines = hits_header + OneOrMore(hitline).setResultsName('hitlines')
        self.hitlines = hitlines

        # parse:
        """
>ref|YP_134121.1| hypothetical protein pNG4041 [Haloarcula marismortui ATCC 43049]
 gb|AAV44415.1| unknown [Haloarcula marismortui ATCC 43049]
          Length = 730

 Score = 27.3 bits (59), Expect = 5.1
 Identities = 12/26 (46%), Positives = 13/26 (50%)
 Frame = -1

Query: 100 CRSCPGQWHQSPEPQAPAVPVQCSRH 23
           CR     WH      APAVP +C RH
Sbjct: 168 CRKAEAGWHW-----APAVPPKCQRH 188


        """

        sequence_name = SkipTo(LineEnd() + "Length =", include=True)
        sequence_name = sequence_name.setResultsName("sequence_name")
        
        alignment_header = Literal(">") + sequence_name + named_int('length')
        self.alignment_header = alignment_header.setResultsName('header')
        
        # parse:
        "39/167 (23%)"
        score_triple = named_int("x") + "/" + named_int("y") + "(" + \
                       named_int("percent") + Word("%),")
        score_triple.setParseAction(make_score_triple)
        
        identities = score_triple.setResultsName('identity')
        positives = score_triple.setResultsName('positives')
        gaps = score_triple.setResultsName('gaps')

        self.score = Literal("Score =") + named_float('bits') + \
                "bits (" + named_int('bits_max') + ")," + \
                Word("Expect()" + nums) + "=" + named_float('expect') + \
                SkipTo("Identities").suppress() + \
                "Identities =" + identities + \
                Optional("Positives =" + positives) + \
                Optional("Gaps =" + gaps) + \
                Optional(("Frame =" + named_frame('frame1') + \
                  Optional("/" + named_frame('frame2'))) |
                         ("Strand =" + strand('frame1') + '/' + strand('frame2')))
                 
        alignment_triple = Literal("Query:").suppress() + \
                           Word(nums).setParseAction(make_int) + \
                           Word(gapped_sequence) + \
                           Word(nums).setParseAction(make_int) + \
                           LineEnd().suppress() + \
                           SkipTo(LineEnd().suppress()) + \
                           Literal("Sbjct:").suppress() + \
                           Word(nums).setParseAction(make_int) + \
                           Word(gapped_sequence) + \
                           Word(nums).setParseAction(make_int)

        chunks = Group(OneOrMore(alignment_triple)).setResultsName('chunks')
        chunks.setParseAction(process_alignment_triple)

        score = Group(self.score + chunks)
        scores = OneOrMore(score).setResultsName('scores')
        
        self.alignments = self.alignment_header + scores
        self.alignments = self.alignments.setResultsName('alignments')
        self.alignments.setParseAction(construct_blast_alignment)

        # parse: (self.query)
        """
Query= EKJFPUJ01A00BV length=100 xy=0304_1625 region=1
run=R_2007_02_22_16_50_08_
         (100 letters)



                                                                 Score    E
Sequences producing significant alignments:                      (bits) Value

ref|YP_136807.1| hypothetical protein rrnAC2273 [Haloarcula mari...    36   0.011
ref|YP_134259.1| hypothetical protein pNG5141 [Haloarcula marism...    28   3.0  
ref|YP_134121.1| hypothetical protein pNG4041 [Haloarcula marism...    27   5.1  
gb|ABI23033.1| succinyl-diaminopimelate desuccinylase [unculture...    27   8.7  

>ref|YP_136807.1| hypothetical protein rrnAC2273 [Haloarcula marismortui ATCC 43049]
 gb|AAV47101.1| unknown [Haloarcula marismortui ATCC 43049]
          Length = 586

 Score = 36.2 bits (82), Expect = 0.011
 Identities = 17/31 (54%), Positives = 20/31 (64%)
 Frame = +1

Query: 4   NNRETSDAYCIELARRGLVVLEIDAIGRGNS 96
           N++ET   + IE ARRG VVL ID  G G S
Sbjct: 77  NSKETQSPFAIEYARRGFVVLAIDQTGHGYS 107


>ref|YP_134259.1| hypothetical protein pNG5141 [Haloarcula marismortui ATCC 43049]
 gb|AAV44553.1| unknown [Haloarcula marismortui ATCC 43049]
          Length = 291

 Score = 28.1 bits (61), Expect = 3.0
 Identities = 12/26 (46%), Positives = 14/26 (53%)
 Frame = -1

Query: 100 CRSCPGQWHQSPEPQAPAVPVQCSRH 23
           CR    +WH      APAVP +C RH
Sbjct: 42  CRETEARWHW-----APAVPPECRRH 62


>ref|YP_134121.1| hypothetical protein pNG4041 [Haloarcula marismortui ATCC 43049]
 gb|AAV44415.1| unknown [Haloarcula marismortui ATCC 43049]
          Length = 730

 Score = 27.3 bits (59), Expect = 5.1
 Identities = 12/26 (46%), Positives = 13/26 (50%)
 Frame = -1

Query: 100 CRSCPGQWHQSPEPQAPAVPVQCSRH 23
           CR     WH      APAVP +C RH
Sbjct: 168 CRKAEAGWHW-----APAVPPKCQRH 188


>gb|ABI23033.1| succinyl-diaminopimelate desuccinylase [uncultured euryarchaeote
           ARMAN-2]
          Length = 291

 Score = 26.6 bits (57), Expect = 8.7
 Identities = 11/24 (45%), Positives = 17/24 (70%)
 Frame = +1

Query: 16  TSDAYCIELARRGLVVLEIDAIGR 87
           TSD   IE+A +G++ L I A+G+
Sbjct: 188 TSDGMSIEIAEKGVLWLRITAVGK 211


"""

        # full query (from start of hits on to next 'Query=')
        self.the_hits = OneOrMore(self.alignments).setResultsName('hits')

        self.query = Group(self.no_hits | (self.hitlines + self.the_hits))
        self.query = self.query.setResultsName('query')

        # full blast output
        self.blast_record = Optional(self.blast_top) + \
                            ((self.query_line + self.db_info) |
                             (Optional(self.db_info) + self.query_line)) + \
                             self.query

        def construct_blast_query_object_with_self(x):
            return construct_blast_query_object(x, parser_obj=self)
        self.blast_record.setParseAction(construct_blast_query_object_with_self)

        self.blast_output = OneOrMore(self.blast_record) + StringEnd()

### Helper/utility data & functions ###

e_val = nums + "e-."
frame = "123-+"
gapped_sequence = alphas + '-*'

def make_comma_int(s):
    return int(s[0].replace(",", ""))

make_int = lambda x: int(x[0])
def make_float(x):
    x = x[0]
    if x.startswith('e'):
        x = '1' + x
    return float(x)

named_int = lambda name :  \
         Word(nums).setParseAction(make_int).copy().setResultsName(name)
named_float = lambda name :  \
         Word(e_val).setParseAction(make_float).copy().setResultsName(name)
named_comma_int = lambda name :  \
         Word(nums+',').setParseAction(make_comma_int).copy().setResultsName(name)
named_frame = lambda name :  \
         Word(frame).setParseAction(make_int).copy().setResultsName(name)

def parse_frame(x):
    if x[0] == 'Plus':
        return 1
    elif x[0] == 'Minus':
        return -1
    assert 0
    
strand = lambda name : \
         (Literal("Plus") | Literal("Minus")).setParseAction(parse_frame).copy().setResultsName(name)
    
def make_score_triple(toks):
    return (ScoreTriple(toks.x, toks.y, toks.percent),)

class ScoreTriple:
    def __init__(self, x, y, percent):
        self.x = x
        self.y = y
        self.percent = percent

    def __repr__(self):
        return 'ScoreTriple<%d, %d, %d>' % (self.x, self.y, self.percent,)

def build_short_sequence_name(name, max_len=20):
    if len(name) < max_len:
        return name

    name_l = name.split()
    if len(name_l) > 1:
        return build_short_sequence_name(name_l[0], max_len)

    name = name_l[0]
    if len(name) > max_len:
        name = name[:max_len-3] + '...'
    return name

def process_hitline(toks):
    toks = toks[0]
            
    name = toks[0]
    return BlastHitSummary(name, toks.total_score, toks.total_expect,
                           toks.N)


def process_alignment_triple(x):
    """
    Package the alignment information into a contiguous triple of strings.
    """
    x = x[0]

    query_seq = subject_seq = alignment = ""
    for i in range(0, len(x), 7):
        (_, q, _, a, _, s, _) = x[i:i+7]
        query_seq += q
        alignment += a
        subject_seq += s
        
    query_start = x[0]
    query_end = x[-5]

    subject_start = x[4]
    subject_end = x[-1]

    return [(query_start, query_end, query_seq,
            subject_start, subject_end, subject_seq,
            alignment)]

def construct_blast_alignment(p):
    matches = []
    for score in p.scores:
        alignment_info = score.chunks
        
        (query_start, query_end, query_seq,
         subject_start, subject_end, subject_seq,
         alignment) = alignment_info[0]
        
        m = BlastSubjectSubmatch(score.expect, score.frame1,
                                 score.frame2, score.bits,
                                 score.bits_max, query_start,
                                 query_end, query_seq, subject_start,
                                 subject_end, subject_seq, alignment)
        
        matches.append(m)

    subject_name = p.header.sequence_name[0]
    subject_name = subject_name.replace("\n", " ")
    
    hits = BlastSubjectHits(subject_name, p.header.length, matches)

    return [hits]

def construct_blast_query_object(x, parser_obj=None):
    """
    Construct the BlastQuery object, containing all information about a
    specific query & results.
    """
    if x.blast_type:
        program = x.blast_type
        version = x.blast_version
        date = x.blast_date
    else:
        program = parser_obj.blast_program
        version = parser_obj.blast_version
        date = parser_obj.blast_date
        

    query_sequence_name = x.query_header.query_name[0]
    query_sequence_name = query_sequence_name.replace("\n", " ")

    query_sequence_length = x.query_header.query_size

    # in some BLAST output files, the database info is only contained
    # in the header, and each BLAST record does *not* have it.  If that's
    # the case (i.e. no x.database_name) retrieve it from the parser object.

    if x.database_name:
        dbname = x.database_name[0].strip()
        database = BlastDatabase(dbname, x.database_num_seqs,
                                 x.database_total_seq)
    else:
        database = parser_obj.blast_database

    hitlines = x.query.hitlines
    hits = x.query.hits

    if hits:
        for (summary, hit) in zip(hitlines, hits):
            hit.total_expect = summary.total_expect
            hit.total_score = summary.total_score
    else:
        hitlines = []
        hits = []

    q = BlastQuery(program, version, date,
                   query_sequence_name, query_sequence_length, database,
                   hitlines, hits)

    return q

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    p = Popen(cmd.split(' '), stdout=PIPE)
    return p.communicate()

def doWork( args ):
    """ Main wrapper"""
    
    for record in parse_file(args.blast_file):
      print '-', record.query_name, record.database.name
      for hit in record.hits:
         print '--', hit.subject_name, hit.subject_length
         print '  ', hit.total_score, hit.total_expect
         for submatch in hit:
            print submatch.expect, submatch.bits
            
            print submatch.query_sequence
            print submatch.alignment
            print submatch.subject_sequence
                
            
            
            
            
            
    
    

    """
# run somethign external in threads
pool = Pool(6)
cmds = ['ls -l', 'ls -alh', 'ps -ef']
print pool.map(runCommand, cmds)
"""

    """
# parse a file
try:
with open(filename, "r") as fh:
for line in fh:
print line
except:
print "Error opening file:", filename, exc_info()[0]
raise
"""

    """
fig = plt.figure()

#-----
# make a 3d plot
ax = fig.add_subplot(111, projection='3d')
ax.scatter(points[:,0],
points[:,1],
points[:,2],
#edgecolors='none',
#c=colors,
#s=2,
#marker='.'
)

#-----
# make a 2d plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(points[:,0],
points[:,1],
'*g')

#-----
# show figure
plt.show()
# or save figure
plt.savefig(filename,dpi=300,format='png')

#-----
# clean up!
plt.close(fig)
del fig
"""

    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--blast_file', help="blast output file")
    #parser.add_argument('input_file2', help="gut_img_ids")
    #parser.add_argument('input_file3', help="oral_img_ids")
    #parser.add_argument('input_file4', help="ids_present_gut_and_oral.csv")
    #parser.add_argument('output_file', help="output file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
