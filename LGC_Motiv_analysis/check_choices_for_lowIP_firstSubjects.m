% extract choices options
choicesOptions = [summary.choiceOptions.monetary_amount.left; summary.choiceOptions.monetary_amount.right];
choicesOptions=choicesOptions(:,strcmp(summary.choiceOptions.R_or_P,'R'));
% extract choice made
choices = summary.choice;
choices = choices(strcmp(summary.choiceOptions.R_or_P,'R'));

nRtrials = 24;
nChoiceOfInterest = 0;
choiceHighMoney = 0;
choiceLowMoney = 0;
for iTrial = 1:nRtrials
    if ismember(choices(iTrial),[-2,-1,1,2])
        if (choicesOptions(1,iTrial) == 1.8 && choicesOptions(2,iTrial) == 1.75) ||...
                (choicesOptions(1,iTrial) == 1.75 && choicesOptions(2,iTrial) == 1.8)
            nChoiceOfInterest =  nChoiceOfInterest + 1;
            if choicesOptions(1,iTrial) == 1.8 && choicesOptions(2,iTrial) == 1.75
                if ismember(choices(iTrial),[-2, -1]) % choice left
                    choiceHighMoney = choiceHighMoney + 1;
                elseif ismember(choices(iTrial),[1,2]) % choice right
                    choiceLowMoney = choiceLowMoney + 1;
                end
            elseif choicesOptions(1,iTrial) == 1.75 && choicesOptions(2,iTrial) == 1.8
                if ismember(choices(iTrial),[-2, -1]) % choice left
                    choiceLowMoney = choiceLowMoney + 1;
                elseif ismember(choices(iTrial),[1,2]) % choice right
                    choiceHighMoney = choiceHighMoney + 1;
                end
            end
        end
    end

end % for loop