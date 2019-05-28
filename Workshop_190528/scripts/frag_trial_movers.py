from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.core.fragment import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.teaching import PyMOLMover
from pyrosetta.rosetta.protocols.moves import *
import numpy as np
from sys import argv
pyrosetta.init(options='-ex1 -ex2aro -mute core.io -mute core.scoring -mute core.pack')


def make_scorefxn():
    # Set up score functions (for full atom and centroid)

    scorefxn = get_fa_scorefxn()
    return scorefxn


def make_Fragmovers(movemap):
    fragset3 = ConstantLengthFragSet(3)
    fragset3.read_fragment_file("aat000_03_05.200_v1_3")
    mover_3mer = ClassicFragmentMover(fragset3, movemap)
    fragset9 = ConstantLengthFragSet(9)
    fragset9.read_fragment_file("aat000_09_05.200_v1_3")
    mover_9mer = ClassicFragmentMover(fragset9, movemap)

    return mover_3mer, mover_9mer


def make_repeatMover(mover_in, n):
    repeat_mover = RepeatMover(mover_in, n)

    return repeat_mover


def make_trial_mover(mover_in, mc):
    trial_mover = TrialMover(mover_in,mc)

    return trial_mover


def main_loop(pose, num_trys):

    score = make_scorefxn()
    print("inital score is " + str(score(pose)) + "\n")
    kt = 1
    mc = MonteCarlo(pose,score,kt)
    movemap = MoveMap()
    movemap.set_bb(True)

    frag_mover3, frag_mover9 = make_Fragmovers(movemap)
    repeat_mover3 = make_repeatMover(frag_mover3,5)
    repeat_mover9 = make_repeatMover(frag_mover9,5)
    trial_mover3 = make_trial_mover(repeat_mover3,mc)
    trial_mover9 = make_trial_mover(repeat_mover9, mc)
    seq_mover = SequenceMover()
    seq_mover.add_mover(trial_mover9)
    seq_mover.add_mover(trial_mover3)

    for i in range(0, num_trys):
        print('attempt #' + str(i) + "\n")
        seq_mover.apply(pose)

    print('final score is ' + str(score(pose)))
    return pose

script, sequence, num_trys = argv

with open(sequence,'r') as f:
    seq = f.read().strip()
    print(seq)

pose = pose_from_sequence(seq, 'fa_standard',auto_termini = True)
pose.dump_pdb("initial_model.pdb")

pmm = PyMOLMover()
pmm.apply(pose)
pmm.keep_history(True)
observer = AddPyMOLObserver(pose, True)

pose = main_loop(pose, int(num_trys))
pose.dump_pdb('final_model.pdb')




