###
## helper functions related to wrangling the covariates, such as conversion to
## dummy matrix format
###

require(caret)
require(tidyverse)

prep_model_data <- function(rds_path) {
    df <- readRDS(rds_path) %>%
        select(-c(education_level, ownership, num_cats)) %>% # waste of time
        mutate_at(
            vars(village:land_for_agriculture, sign_animals:windows_in_bedroom),
            ~as_factor(str_replace_all(.x, ' ', '_'))
        )
    
    levels(df$house_age) <- c('2_to_6', '7_or_more', 'one_or_less')
    df$truth <- df$infestation
    return(df)
}

dummify_vars <- function(df_full) {
    disc <- df_full %>%
        select(-c(id:village) & where(is.factor)) %>%
        modify_if(~any(levels(.x) == 'Other'), fct_relevel, 'Other')
    
    dv <- dummyVars(
        ~ .,
        disc,
        sep = '.',
        fullRank = TRUE
    )
    ## only return vars nec. for kriging:
    bind_cols(
        select(df_full, long, lat),
        as_tibble(predict(dv, newdata = disc), .name_repair = 'minimal'),
        select(df_full, contains('num_'))
    )
}

undummify <- function(df_ltfr) {
    df_ltfr %>%
        pivot_longer(contains('.')) %>%
        separate(name, into = c('f', 'l'), sep = '\\.') %>%
        filter(value == 1) %>%
        pivot_wider(!c(f, l, value), names_from = f, values_from = l) %>%
        relocate(infestation, truth, .after = last_col()) # need at end for make_formula
}

rescale_cont_vars <- function(df) {
    df_cont <- select(df, contains('num_'), dist_perim, density)
    df_oth <- select(df, !c(contains('num_'), dist_perim, density))

    bind_cols(
        df_oth,
        predict(
            suppressWarnings(
                preProcess(
                    df_cont,
                    method = c('center', 'scale')
                )
            ),
            newdata = df_cont
        )
    ) %>%
        ## put before inf for make_formula
        relocate(contains('num_'), dist_perim, density, .before = infestation)
}

make_dummy_mat <- function(df_org, as_tibble=FALSE, contr='contr.bayes',
                           rescale = TRUE, inter = TRUE, trim = FALSE) {
    ## wrapper to get cov matrix with consistent reference levels
    ## if as_tibble true, will keep geospatial identifiers and infestation status
    require(bayestestR)

    ## make sure contrasts are correctly set
    options(contrasts = c(contr, 'contr.poly'))

    df_trim <- df_org[ ,-c(1:4)]
    if (trim) {
        df_trim <- df_trim %>%
            select_if(~ length(unique(.x)) > 1)
    }
    
    dat_disc <- df_trim %>%
        select_if(is.factor) %>%
        mutate(infestation = df_org$infestation)
    
    dat_cont <- df_trim %>%
        select(-infestation) %>%
        select_if(is.double)

    if (rescale) {
        dat_cont <- predict(
            preProcess(
                dat_cont,
                method = c('center', 'scale')
            ),
            newdata = dat_cont
        )
    }
    
    ## relevel/collapse each factor to remove any low level counts
    dat_disc$bed_hygiene <- dat_disc$bed_hygiene %>%
        relevel('buena')
    
    dat_disc$bird_nests_inside <- dat_disc$bird_nests_inside %>%
        relevel('no')
    
    dat_disc$chicken_coop_location <- dat_disc$chicken_coop_location %>%
        relevel('no')
    
    dat_disc$bedroom_clutter <- dat_disc$bedroom_clutter %>%
        relevel('no acumula')
    
    dat_disc$dark_main_room <- dat_disc$dark_main_room %>%
        relevel('no')
    
    dat_disc$education_level <- dat_disc$education_level %>%
        relevel('ninguno')
    
    dat_disc$floor_type <- dat_disc$floor_type %>%
        relevel('ladrillo')
    
    dat_disc$construction_pile_type <- dat_disc$construction_pile_type %>%
        relevel('ninguno')
    
    dat_disc$firewood_location <- dat_disc$firewood_location %>%
        relevel('ninguno')
    
    dat_disc$grain_storage_in_house <- dat_disc$grain_storage_in_house %>%
        relevel('no')

    if (any(levels(dat_disc$house_age) == 'menos de 1 año')) {
        dat_disc$house_age <- dat_disc$house_age %>%
            relevel('menos de 1 año')
    }
    
    dat_disc$house_hygiene <- dat_disc$house_hygiene %>%
        relevel('buena')
    
    dat_disc$kitchen_location <- dat_disc$kitchen_location %>%
        relevel('shared_none')
    
    dat_disc$land_for_agriculture <- dat_disc$land_for_agriculture %>%
        relevel('no')
    
    dat_disc$material_house_wall <- dat_disc$material_house_wall %>%
        relevel('brick_block_oth')
    
    dat_disc$material_roof <- dat_disc$material_roof %>%
        relevel('alum_cement')
    
    dat_disc$ownership <- dat_disc$ownership %>%
        relevel('owned')
    
    dat_disc$sign_animals <- dat_disc$sign_animals %>%
        relevel('FALSE')
    
    dat_disc$sign_rats <- dat_disc$sign_rats %>%
        relevel('FALSE')
    
    dat_disc$windows_in_bedroom <- dat_disc$windows_in_bedroom %>%
        relevel('FALSE')
    
    dat_disc$condition_house_wall <- dat_disc$condition_house_wall %>%
        relevel('buen estado')
    
    dat_disc$condition_bedroom_wall <- dat_disc$condition_bedroom_wall %>%
        relevel('buen estado')

    m_orth <- model.matrix(infestation ~ . - 1, data=dat_disc)[,-1]
    if (inter)
        cov_matrix <- cbind(inter = 1, m_orth, as.matrix(dat_cont))
    else
        cov_matrix <- cbind(m_orth, as.matrix(dat_cont))
    
    if (as_tibble) {
        cov_matrix %<>%
            as_tibble(.name_repair='minimal') %>%
            bind_cols(df_org[,1:4], ., infestation=df_org$infestation)
    }
    
    return(cov_matrix)
}
