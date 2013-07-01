oxygen
======

Proton-proton correlations in oxygen-proton interactions






- LednicCorrel
        Makefile -> lednic.exe
        v root-e  auto_lednic.C -> led_ZR=X.XX

- LednicCorrel/MoreInput
        Makefile -> lednic_input.exe
        v root-e auto_lednic_input.C -> rozne obrazky s roznymi parametrami

- Experiment/DST
        Makefile -> redact.exe
        redact.exe -> oxygen_DST.asc
        v root-e oxygen_DST.C -> oxygen_DST.root (link na tento subor)

- Model/Moscow
        v root-e moscow_model.C -> moscow_model.root (link na tento subor)

- Model/Musul
        Makefile -> musul_model.exe
        musul_model.exe -> musul_model.asc
        v root-e musul_model.C -> musul_model.root (link na tento subor)

- Model/ModelCorrect
        !!! vsetko len pre protony v impulsnom intervale (1.75, 4.75) !!!
        .x find_error_landau.C -> obrazok s exp. chybami protonov delta_p/tot_p
        .x correct_on_landau.C(1) -> moscow_model_c.root (link na tento subor)
        .x correct_on_landau.C(2) -> musul_model_c.root (link na tento subor)

- Model
        .x correlation_model.C -> korelacna funkcia pre model, model poopraveny
        na chybu correct_on_moscow.asc alebo correct_on_musul.asc

- Experiment
        .x correlation_DST.C
