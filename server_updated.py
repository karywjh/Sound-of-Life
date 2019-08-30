from flask import Flask, request, render_template

# Initialize Flask app
app = Flask(__name__, static_folder="static")

import os
import urllib.request
from app import app
from flask import Flask, flash, request, redirect, render_template
from werkzeug.utils import secure_filename
from EasyMIDI import EasyMIDI,Track,Note,Chord,RomanChord
from random import choice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from midi2audio import FluidSynth
import time

ALLOWED_EXTENSIONS = set(['txt','fasta',"fst"])

def rco(seq) :
    #coding_dna = Seq('', IUPAC.unambiguous_dna)
    coding_dna = Seq(seq, IUPAC.unambiguous_dna)
    basepair = str(coding_dna)
    template_dna = coding_dna.reverse_complement()
    mrna = coding_dna.transcribe()
    aminoacids = mrna.translate(table = 2, to_stop=True)
    aminomusicsequence = aminoacids
    print(basepair)

    duratlist = []
    for i in range(len(basepair)):
        if basepair[i] == 'A':
            duratlist.append(1/2)
        elif basepair[i] == 'T':
            duratlist.append(1/4)
        elif basepair[i] == 'C':
            duratlist.append(1/8)
        elif basepair[i] == 'G':
            duratlist.append(3/8)

    output = []
    for i, val in enumerate(duratlist):
        if (i + 1) % 3 == 0:
            output.append(1 - (duratlist[i-1] + duratlist[i-2]))
        else:
            output.append(val)

    chunks = []
    temp_chunk = []
    for val in output:
        temp_chunk.append(val)
        if len(temp_chunk) == 3:
            chunks.append(temp_chunk)
            temp_chunk = []
    return chunks

def generate(dna):
    DNA=Seq(dna,IUPAC.unambiguous_dna)
    seq=DNA.translate()
    seq=str(seq)

    durat_list = rco(dna)
    print(durat_list)

    # Generation

    easyMIDI = EasyMIDI(numTracks=3)
    track1 = Track('acoustic grand piano')
    #track2 = Track('FX 5 (brightness)')
    track2 = Track('Whistle')
    track3 = Track('String Ensemble 1')

    key_input = 'C'

    c2 = Note('C', octave = 2, duration = 1/8, volume = 100)
    d2 = Note('D', 2, 1/8)
    e2 = Note('E', 2, 1/8)
    eflat2 = Note('Eb', 2, 1/8)
    f2 = Note('F', 2, 1/8)
    fsharp2 = Note('F#', 2, 1/8)
    g2 = Note('G', 2, 1/8)
    a2 = Note('A', 2, 1/8)
    b2 = Note('B', 2, 1/8) 
    bflat2 = Note('Bb',2,1/8)

    #3rd Octave
    c = Note('C', octave = 3, duration = 1/8, volume = 100)
    d = Note('D', 3, 1/8)
    e = Note('E', 3, 1/8)
    eflat = Note('Eb',3, 1/8)
    f = Note('F', 3, 1/8)
    fsharp = Note('F#', 3, 1/8)
    g = Note('G', 3, 1/8)
    a = Note('A', 3, 1/8)
    b = Note('B', 3, 1/8)
    bflat = Note('Bb',3,1/8)

    #duration=3/8, 3rd Octave
    d_38 = Note('D', 3, 3/8)
    fsharp_38 = Note('F#',3, 3/8)
    a_38 = Note('A', 3, 3/8)
    c4_38 = Note('C', 3, 3/8)

    #duration=1/4, 3rd Octave
    d_14 = Note('D', 3, 1/4)
    fsharp_14 = Note('F#', 3, 1/4)
    a_14 = Note('A', 3, 1/4)
    c4_14 = Note('C', 3, 1/4)

    #4th Octave
    c4 = Note('C', octave = 4, duration = 1/8, volume = 100)
    d4 = Note('D', 4,1/8)
    e4 = Note('E', 4,1/8)
    eflat4 = Note('Eb',4, 1/8)
    f4 = Note('F', 4,1/8)
    fsharp4 = Note('F#',4,1/8)
    g4 = Note('G', 4,1/8)
    a4 = Note('A', 4,1/8)
    b4 = Note('B', 4,1/8)
    bflat4 = Note('Bb',4,1/8)

    #5th Octave
    c5 = Note('C',5,1/8)
    d5 = Note('D',5,1/8)
    e5 = Note('E',5,1/8)
    eflat5 = Note('Eb',5, 1/8)
    f5 = Note('F', 5,1/8)
    fsharp5 = Note('F#',5,1/8)
    g5 = Note('G', 5,1/8)
    a5 = Note('A', 5,1/8)
    b5 = Note('B', 5,1/8)
    bflat5 = Note('Bb',5,1/8)

    #6th Octave
    c6 = Note('C',6,1/8)
    d6 = Note('D',6,1/8)
    e6 = Note('E',6,1/8)
    eflat6 = Note('Eb',6, 1/8)
    f6 = Note('F', 6,1/8)
    fsharp6 = Note('F#',6,1/8)
    g6 = Note('G', 6,1/8)
    a6 = Note('A', 6,1/8)
    b6 = Note('B', 6,1/8)
    bflat6 = Note('Bb',6,1/8)

    #chordC = Chord([c,e,g])  #-- a chord of notes C, E and G

    # Track3
    # abcbcbcb pattern
    Idiv = [c,e,g,e,g,e,g,e]
    I6div = [e,g,c4,g,c4,g,c4,g]
    I64div = [g2,c,e,c,e,c,e,c]
    Idom7div = [c,g,bflat,g,bflat,g,bflat,g]

    IIdiv = [d,f,a,f,a,f,a,f]
    II6div = [f2,a2,d,a2,d,a2,d,a2]
    II64div = [a2,d,f,d,f,d,f,d]

    IIIdiv = [e,g,b,g,b,g,b,g]

    IVdiv = [f,a,c4,a,c4,a,c4,a]
    IV6div = [a2,c,f,c,f,c,f,c]
    IV64div = [c,f,a,f,a,f,a,f]
    IV7div = [f2,c,eflat,c,eflat,c,eflat,c]

    Vdiv = [g2,b2,d,b2,d,b2,d,b2]
    V6div = [b2,d,g,d,g,d,g,d]
    V64div = [c,f,a,f,a,f,a,f]
    V7div = [g2,d,f,d,f,d,f,d]
    VofVdiv = [d,fsharp,c4,fsharp,c4,fsharp,c4,fsharp]

    VIdiv = [a2,c,e,c,e,c,e,c]
    VI6div = [c,e,a,e,a,e,a,e]

    def set_Duration(chordname, octave, durat):
        cd = RomanChord(chordname, octave, durat, key_input, True, 100)
        return cd

    #regular roman chords
    I = RomanChord(numeral='I', octave=4, duration=0.25, key=key_input, major=True, volume=100)
    II = RomanChord (numeral='II', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    III = RomanChord (numeral='III', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    IV = RomanChord(numeral='IV', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    V = RomanChord(numeral='V', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    VI = RomanChord(numeral='VI', octave=3, duration=0.25, key=key_input, major=True, volume=100)

    I6 = RomanChord(numeral='I*', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    I64 = RomanChord(numeral='I**', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    Idom7 = RomanChord(numeral='Idom7', octave=3, duration=0.25, key=key_input, major=True, volume=100)

    II6 = RomanChord(numeral='II*', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    II64 = RomanChord(numeral='II**', octave=3, duration=0.25, key=key_input, major=True, volume=100)

    IV6 = RomanChord(numeral='IV*', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    IV64 = RomanChord(numeral='IV**', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    IV7 = RomanChord(numeral='IVdom7', octave=3, duration=0.25, key=key_input, major=True, volume=100)

    V6 = RomanChord(numeral='V*', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    V64 = RomanChord(numeral='V**', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    V7 = RomanChord(numeral='V7', octave=3, duration=0.25, key=key_input, major=True, volume=100)
    # VofV = Chord([d,fsharp,a,c4])
    # VofV_38 = Chord([d_38,fsharp_38,a_38,c4_38])
    # VofV_14 = Chord([d_14,fsharp_14,a_14,c4_14])

    VI6 = RomanChord(numeral='VI*', octave=3, duration=0.25, key=key_input, major=True, volume=100)

    PAC_cadence = [V6,I]

    chordPredict =[]

    def Starting_Chord (track):
        track.addChord(set_Duration('I',3,3/8))
        track.addChord(set_Duration('I',3,1/8))
        track.addChord(set_Duration('I',3,1/4))
        track.addChord(set_Duration('I',3,1/4))

    for letter in seq:
        if letter=='S':
            #I
            chordPredict.append('I')
            track2.addChord(set_Duration('I',4,3/8))
            track2.addChord(set_Duration('I',4,1/8))
            track2.addChord(set_Duration('I',4,1/4))
            track2.addChord(set_Duration('I',4,1/4))
            track3.addNote(Idiv)
        elif letter=='L':
            #V
            chordPredict.append('V')
            track2.addChord(set_Duration('V',3,3/8))
            track2.addChord(set_Duration('V',3,1/8))
            track2.addChord(set_Duration('V',3,1/4))
            track2.addChord(set_Duration('V',3,1/4))
            track3.addNote(Vdiv)
        elif letter=='A':
            #IV
            chordPredict.append('IV')
            track2.addChord(set_Duration('IV',3,3/8))
            track2.addChord(set_Duration('IV',3,1/8))
            track2.addChord(set_Duration('IV',3,1/4))
            track2.addChord(set_Duration('IV',3,1/4))
            track3.addNote(IVdiv)
        elif letter=='G':
            #V7
            chordPredict.append('V7')
            track2.addChord(set_Duration('V7',3,3/8))
            track2.addChord(set_Duration('V7',3,1/8))
            track2.addChord(set_Duration('V7',3,1/4))
            track2.addChord(set_Duration('V7',3,1/4))
            track3.addNote(V7div)
        elif letter=='K':
            #VI
            chordPredict.append('VI')
            track2.addChord(set_Duration('VI',3,3/8))
            track2.addChord(set_Duration('VI',3,1/8))
            track2.addChord(set_Duration('VI',3,1/4))
            track2.addChord(set_Duration('VI',3,1/4))
            track3.addNote(VIdiv)
        elif letter=='T':
            #I6
            chordPredict.append('I*')
            track2.addChord(set_Duration('I*',3,3/8))
            track2.addChord(set_Duration('I*',3,1/8))
            track2.addChord(set_Duration('I*',3,1/4))
            track2.addChord(set_Duration('I*',3,1/4))
            track3.addNote(I6div)
        elif letter=='D':
            #V6
            chordPredict.append('V*')
            track2.addChord(set_Duration('V*',3,3/8))
            track2.addChord(set_Duration('V*',3,1/8))
            track2.addChord(set_Duration('V*',3,1/4))
            track2.addChord(set_Duration('V*',3,1/4))
            track3.addNote(V6div)
        elif letter=='E':
            #IV6
            chordPredict.append('IV*')
            track2.addChord(set_Duration('IV*',3,3/8))
            track2.addChord(set_Duration('IV*',3,1/8))
            track2.addChord(set_Duration('IV*',3,1/4))
            track2.addChord(set_Duration('IV*',3,1/4))
            track3.addNote(IV6div)
        elif letter=='P':
            #Idom7
            chordPredict.append('Idom7')
            track2.addChord(set_Duration('Idom7',3,3/8))
            track2.addChord(set_Duration('Idom7',3,1/8))
            track2.addChord(set_Duration('Idom7',3,1/4))
            track2.addChord(set_Duration('Idom7',3,1/4))
            track3.addNote(Idom7div)
        elif letter=='N':
            #VI6
            chordPredict.append('VI*')
            track2.addChord(set_Duration('VI*',3,3/8))
            track2.addChord(set_Duration('VI*',3,1/8))
            track2.addChord(set_Duration('VI*',3,1/4))
            track2.addChord(set_Duration('VI*',3,1/4))
            track3.addNote(VI6div)
        elif letter=='R':
            #II6
            chordPredict.append('II*')
            track2.addChord(set_Duration('II*',3,3/8))
            track2.addChord(set_Duration('II*',3,1/8))
            track2.addChord(set_Duration('II*',3,1/4))
            track2.addChord(set_Duration('II*',3,1/4))
            track3.addNote(II6div)
        elif letter=='F':
            #VofV
            #VofV = Chord([d,fsharp,a,c4])
            # chordPredict.append('V/V')
            # track2.addChord(VofV_38)
            # track2.addChord(VofV)
            # track2.addChord(VofV_14)
            # track2.addChord(VofV_14)
            # track3.addNote(VofVdiv)
            #V
            chordPredict.append('V')
            track2.addChord(set_Duration('V',3,3/8))
            track2.addChord(set_Duration('V',3,1/8))
            track2.addChord(set_Duration('V',3,1/4))
            track2.addChord(set_Duration('V',3,1/4))
            track3.addNote(Vdiv)
        elif letter=='I':
            #I64
            chordPredict.append('I**')
            track2.addChord(set_Duration('I**',3,3/8))
            track2.addChord(set_Duration('I**',3,1/8))
            track2.addChord(set_Duration('I**',3,1/4))
            track2.addChord(set_Duration('I**',3,1/4))
            track3.addNote(I64div)
        elif letter=='Q':
            #V64
            chordPredict.append('V**')
            track2.addChord(set_Duration('V**',3,3/8))
            track2.addChord(set_Duration('V**',3,1/8))
            track2.addChord(set_Duration('V**',3,1/4))
            track2.addChord(set_Duration('V**',3,1/4))
            track3.addNote(V64div)
            #I
            chordPredict.append('I')
            track2.addChord(set_Duration('I',4,3/8))
            track2.addChord(set_Duration('I',4,1/8))
            track2.addChord(set_Duration('I',4,1/4))
            track2.addChord(set_Duration('I',4,1/4))
            track3.addNote(Idiv)
        elif letter=='C':
            #IV7
            chordPredict.append('IVdom7')
            track2.addChord(set_Duration('IVdom7',3,3/8))
            track2.addChord(set_Duration('IVdom7',3,1/8))
            track2.addChord(set_Duration('IVdom7',3,1/4))
            track2.addChord(set_Duration('IVdom7',3,1/4))
            track3.addNote(IV7div)
        elif letter=='Y':
            #IV64
            chordPredict.append('IV**')
            track2.addChord(set_Duration('IV**',3,3/8))
            track2.addChord(set_Duration('IV**',3,1/8))
            track2.addChord(set_Duration('IV**',3,1/4))
            track2.addChord(set_Duration('IV**',3,1/4))
            track3.addNote(IV64div)
        elif letter=='H':
            #I
            chordPredict.append('I')
            track2.addChord(set_Duration('I',4,3/8))
            track2.addChord(set_Duration('I',4,1/8))
            track2.addChord(set_Duration('I',4,1/4))
            track2.addChord(set_Duration('I',4,1/4))
            track3.addNote(Idiv)
            #V6
            chordPredict.append('V*')
            track2.addChord(set_Duration('V*',3,3/8))
            track2.addChord(set_Duration('V*',3,1/8))
            track2.addChord(set_Duration('V*',3,1/4))
            track2.addChord(set_Duration('V*',3,1/4))
            track3.addNote(V6div)
            #VI
            chordPredict.append('VI')
            track2.addChord(set_Duration('VI',3,3/8))
            track2.addChord(set_Duration('VI',3,1/8))
            track2.addChord(set_Duration('VI',3,1/4))
            track2.addChord(set_Duration('VI',3,1/4))
            track3.addNote(VIdiv)
            #IV6
            chordPredict.append('IV*')
            track2.addChord(set_Duration('IV*',3,3/8))
            track2.addChord(set_Duration('IV*',3,1/8))
            track2.addChord(set_Duration('IV*',3,1/4))
            track2.addChord(set_Duration('IV*',3,1/4))
            track3.addNote(IV6div)
            #I
            chordPredict.append('I')
            track2.addChord(set_Duration('I',4,3/8))
            track2.addChord(set_Duration('I',4,1/8))
            track2.addChord(set_Duration('I',4,1/4))
            track2.addChord(set_Duration('I',4,1/4))
            track3.addNote(Idiv)
        elif letter=='M':
            #II64
            chordPredict.append('II**')
            track2.addChord(set_Duration('II**',3,3/8))
            track2.addChord(set_Duration('II**',3,1/8))
            track2.addChord(set_Duration('II**',3,1/4))
            track2.addChord(set_Duration('II**',3,1/4))
            track3.addNote(II64div)
        elif letter=='W':
            #II6
            chordPredict.append('II*')
            track2.addChord(set_Duration('II*',3,3/8))
            track2.addChord(set_Duration('II*',3,1/8))
            track2.addChord(set_Duration('II*',3,1/4))
            track2.addChord(set_Duration('II*',3,1/4))
            track3.addNote(II6div)
            #V
            chordPredict.append('V')
            track2.addChord(set_Duration('V',3,3/8))
            track2.addChord(set_Duration('V',3,1/8))
            track2.addChord(set_Duration('V',3,1/4))
            track2.addChord(set_Duration('V',3,1/4))
            track3.addNote(Vdiv)
            #I
            chordPredict.append('I')
            track2.addChord(set_Duration('I',4,3/8))
            track2.addChord(set_Duration('I',4,1/8))
            track2.addChord(set_Duration('I',4,1/4))
            track2.addChord(set_Duration('I',4,1/4))
            track3.addNote(Idiv)
        elif letter =='*':
            continue


    getNamedNote = lambda note: note.name + str(note.octave)
    getManyNamedNotes = lambda notes: list(map(getNamedNote, notes))

    # print(chordPredict[0])
    # print(getManyNamedNotes(RomanChord(chordPredict[0]).getNotes()))

    def planB (track, chordPredict, durat_list):
        for i, val in enumerate(chordPredict):
            if i >= len(durat_list):
                break
            N = getManyNamedNotes(RomanChord(val).getNotes())
            for k, val2 in enumerate(durat_list[i]):
                if val2 ==0:
                    break
                m = Note(N[0][:-1], octave=int(N[0][-1]), duration=val2, volume=100)
                p = Note(N[1][:-1], octave=int(N[1][-1]), duration=val2, volume=100)
                q = Note(N[2][:-1], octave=int(N[2][-1]), duration=val2, volume=100)
                mQ = Chord([m,q])
                if dna[i*3+k] == 'A':
                    track.addNote(m)
                elif dna[i*3+k] == 'T':
                    track.addNote(p)
                elif dna[i*3+k] == 'C':
                    track.addNote(q)
                elif dna[i*3+k] == 'G':
                    track.addChord(mQ)

    # Starting_Chord(track1)
    # Starting_Chord(track2)
    # Starting_Chord(track3)
    track2.addChord(PAC_cadence)
    planB(track1, chordPredict, durat_list)
    track3.addChord(PAC_cadence)
    track1.addChord(PAC_cadence)

    easyMIDI.addTrack(track1)
    easyMIDI.addTrack(track2)
    easyMIDI.addTrack(track3)
    easyMIDI.writeMIDI('output.mid')

    fs = FluidSynth(sound_font='SGM.sf2', sample_rate=22050)
    fs.midi_to_audio('output.mid', 'static/output.wav')

def allowed_file(filename):
	return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS
	
@app.route('/', methods=['GET'])
def upload_form():
	return render_template('upload.html')

@app.route('/', methods=['POST'])
def upload_file():
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash('No file selected for uploading')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            aaaa=open(filepath,"r")
            dna = "".join(aaaa.readlines())
            #flash('File successfully uploaded')
            generate(dna)
            flash(dna)
            return redirect('/')
        else:
            flash('Allowed file types are txt, fasta, fst')
            return redirect(request.url)

    return render_template("upload.html")


if __name__ == '__main__':
    app.run()
