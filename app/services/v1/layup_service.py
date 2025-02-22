from typing import List

def parse_layup_sequence(value: str) -> List[float]:
    layups = []
    base_layup = ""
    r_before = 1
    symmetry = False
    r_after = 1

    if "]" in value and value.split("]")[1] != "":
        base_layup = value.split("]")[0].replace("[", "")
        msn = value.split("]")[1]
        if 's' in msn:
            # symmetry
            symmetry = True
            r = msn.split('s')
            if len(r) == 2:
                if r[0] != "":
                    r_before_temp = int(r[0])
                    r_before = r_before_temp
                if r[1] != "":
                    r_after_temp = int(r[1])
                    r_after = r_after_temp
        else:
            # not symmetry
            symmetry = False
            r_before_temp = int(msn)
            r_before = r_before_temp
    else:
        base_layup = value.replace("[", "").replace("]", '')

    for angle_string in base_layup.split('/'):
        angle = float(angle_string)
        layups.append(angle)

    layups_temp = layups.copy()

    for _ in range(1, r_before):
        for layup in layups_temp:
            layups.append(layup)

    layups_temp = layups.copy()

    if symmetry:
        layups_temp_reversed = list(reversed(layups_temp))
        for layup in layups_temp_reversed:
            layups.append(layup)

    layups_temp = layups.copy()

    for _ in range(1, r_after):
        for layup in layups_temp:
            layups.append(layup)

    return layups