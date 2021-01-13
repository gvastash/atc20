rm bscore; for a in `ls -d bsuite/*`; do ./judge_B ./main < $a >> bscore; done
