from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.core.fragment import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.teaching import PyMOLMover
from pyrosetta.rosetta.protocols.moves import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
import numpy as np
from sys import argv

pyrosetta.init(options='-ex1 -ex2aro -mute core.io -mute core.scoring -mute core.pack')


def make_switchers():
    # Set up the switch movers to move between centroid and full atom

    switch_low = SwitchResidueTypeSetMover("centroid")
    switch = SwitchResidueTypeSetMover("fa_standard")
    return switch_low, switch


def make_scorefxn(type):
    # Set up score functions (for full atom and centroid)
    if type == 1:
        scorefxn = get_fa_scorefxn()
    elif type == 2:
        scorefxn = scorefxn_low = create_score_function('score3')

    return scorefxn


def make_Fragmovers(movemap):
    # Load frag sets and make movers

    fragset3 = ConstantLengthFragSet(3)
    fragset3.read_fragment_file("aat000_03_05.200_v1_3")
    mover_3mer = ClassicFragmentMover(fragset3, movemap)
    fragset9 = ConstantLengthFragSet(9)
    fragset9.read_fragment_file("aat000_09_05.200_v1_3")
    mover_9mer = ClassicFragmentMover(fragset9, movemap)

    return mover_3mer, mover_9mer


def make_small_shear(movemap, kt, nmoves):
    small_mover = SmallMover(movemap, kt, nmoves)
    shear_mover = ShearMover(movemap, kt, nmoves)

    return small_mover, shear_mover


def make_repeatMover(mover_in, n):
    repeat_mover = RepeatMover(mover_in, n)

    return repeat_mover


def make_trial_mover(mover_in, mc):
    trial_mover = TrialMover(mover_in, mc)

    return trial_mover


def main_loop(pose, num_trys, pmm):
    # Set up initialisation parameters (kt, score and monte carlo objects)
    score = make_scorefxn(2)
    score_high = make_scorefxn(1)

    switch_low, switch_high = make_switchers()
    switch_low.apply(pose)


    print("inital score is " + str(score(pose)) + "\n")
    kt = 0.1
    mc = MonteCarlo(pose, score, kt)
    mc2 = MonteCarlo(pose, score, kt)

    # Create movemap for controlling the folding
    movemap = MoveMap()
    movemap.set_bb(True)

    # create all of the necessary movers
    frag_mover3, frag_mover9 = make_Fragmovers(movemap)
    small_mover, shear_mover = make_small_shear(movemap, kt, 1)

    # Set up the complex movers
    smallShear_mover = SequenceMover()
    smallShear_mover.add_mover(small_mover)
    smallShear_mover.add_mover(shear_mover)
    repeat_smallShear = make_repeatMover(smallShear_mover, 100)
    smallShear_trial = make_trial_mover(repeat_smallShear, mc2)

    # Set up the insertion/small shear mover
    insertSS_mover_seq = SequenceMover()
    insertSS_mover_seq.add_mover(frag_mover9)
    insertSS_mover_seq.add_mover(smallShear_trial)
    trial_insertSS = make_trial_mover(insertSS_mover_seq, mc)

    # Set up final trial moving that will actually perform folding

    for i in range(0, num_trys):
        print('attempt #' + str(i) + "\n")
        trial_insertSS.apply(pose)
        if i < 100:
            switch_high.apply(pose)
            pmm.apply(pose)
            switch_low.apply(pose)
        elif i % 100 == 0:
            switch_high.apply(pose)
            pmm.apply(pose)
            switch_low.apply(pose)

    switch_high.apply(pose)
    print('final score is ' + str(score_high(pose)))
    return pose


script, sequence, num_trys = argv

with open(sequence, 'r') as f:
    seq = f.read().strip()
    print(seq)

pose = pose_from_sequence(seq, 'fa_standard', auto_termini=True)
pose.dump_pdb("initial_model.pdb")

pmm = PyMOLMover()
pmm.apply(pose)
pmm.keep_history(True)
#observer = AddPyMOLObserver(pose, True)

pose = main_loop(pose, int(num_trys), pmm)
pose.dump_pdb('final_model.pdb')
