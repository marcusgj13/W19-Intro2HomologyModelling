from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from sys import argv

pyrosetta.init(options='-ex1 -ex2aro -mute core.io -mute core.scoring -mute core.pack')
# Implementation of monte carlo folding using movers


def make_switchers():
    # Set up the switch movers to move between centroid and full atom

    switch_low = SwitchResidueTypeSetMover("centroid")
    switch = SwitchResidueTypeSetMover("fa_standard")
    return switch_low, switch


def make_scorefxn():
    # Set up score functions (for full atom and centroid)

    scorefxn_low = create_score_function('score3')
    scorefxn = get_fa_scorefxn()
    return scorefxn_low, scorefxn


def make_small_shear(movemap, kt, nmoves):
    small_mover = SmallMover(movemap, kt, nmoves)
    shear_mover = ShearMover(movemap, kt, nmoves)

    return small_mover, shear_mover


def make_min_mover(scorefxn, movemap):
    min_mover = MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(scorefxn)

    return min_mover


# run the folding algorithm (can run this a bunch of ways)
def monteCarloFold(Input_P, nIterations, trial):

    P = Pose()
    P.assign(Input_P)

    pmm = PyMOLMover()
    pmm.apply(P)
    pmm.keep_history(False)


    switch_low, switch_high = make_switchers()
    scorefxn_low, scorefxn = make_scorefxn()

    # Switch to centroids for simulation
    switch_low.apply(P)
    initial_score = scorefxn_low(P)

    print("Initial score is %f\n" % initial_score)

    # Set up simulation parameters and MonteCarlo mover
    k_t = 1

    # Set up the monteCarlo checker
    mc = MonteCarlo(P, scorefxn_low, k_t)
    movemap = MoveMap()
    movemap.set_bb(True)

    # Set up movers
    small_mover, shear_mover = make_small_shear(movemap,k_t, 5)
    min_mover = make_min_mover(scorefxn_low, movemap)

    for i in range(1, nIterations):

        small_mover.apply(P)
        min_mover.apply(P)
        shear_mover.apply(P)
        min_mover.apply(P)

        mc.boltzmann(P)
        mc.show_scores()
        pmm.apply(P)

    final_score = scorefxn_low(P)
    switch_high.apply(P)
    print("The final score is %f\n" % final_score)
    P.dump_pdb('LowestEnergyModel' + str(trial) + '.pdb')
    return final_score


script, nIterations, nTrials, sequence = argv

with open(sequence,'r') as f:
    seq = f.read().strip()
    print(seq)

# create intial model

Start_P = pose_from_sequence(seq, 'fa_standard',auto_termini = True)
Start_P.dump_pdb('Initial_model.pdb')
all_scores = []


for trial in range(1,int(nTrials)+1):
    score = monteCarloFold(Start_P, int(nIterations), trial)
    if trial == 1:
        with open('folding_results.txt','w') as f:
            f.write('Attempt ' + '\t')
            f.write('Score ' +  '\n')
            f.write(str(trial) + '\t')
            f.write(str(score) + '\n')
    else:
        with open('folding_results.txt','a') as f:
            f.write(str(trial) + '\t')
            f.write(str(score) + '\n')




