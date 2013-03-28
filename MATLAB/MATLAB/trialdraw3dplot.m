
xlabel({'Larval Stage Survival'})
ylabel({'Egg Stage Survival'})
subplot(3,2,1); contourf (EM, LM,ELAdultBeesF1);title({'AdultBees Number at First November'});xlabel({'Larval Stage Survival'});ylabel({'Egg Stage Survival'})
%subplot(3,2,2); contourf (EM, LM,ELHoneystoreW1); title({'Honeystorage at First End of Winter'});xlabel({'Larval Stage Survival'});ylabel({'Egg Stage Survival'})
subplot(3,2,3); contourf (EM, LM,ELAdultBeesF2); title({'AdultBees Number at Second November'});xlabel({'Larval Stage Survival'});ylabel({'Egg Stage Survival'})
%subplot(3,2,4); contourf (EM, LM,ELHoneystoreW2); title({'Honeystorage at Second End of Winter'});xlabel({'Larval Stage Survival'});ylabel({'Egg Stage Survival'})
subplot(3,2,5); contourf (EM, LM,ELAdultBeesF3); title({'AdultBees Number at Third November'});xlabel({'Larval Stage Survival'});ylabel({'Egg Stage Survival'})
%subplot(3,2,6); contourf (EM, LM,ELHoneystoreW3); title({'Honeystorage at Third End of Winter'});xlabel({'Larval Stage Survival'});ylabel({'Egg Stage Survival'})
subplot(3,2,2); contourf (EM, LM,ELAdultBeesW1); title({'AdultBees Number at First End of Winter'});xlabel({'Larval Stage Survival'});ylabel({'Egg Stage Survival'})
subplot(3,2,4); contourf (EM, LM,ELAdultBeesW2); title({'AdultBees Number at Second End of Winter'});xlabel({'Larval Stage Survival'});ylabel({'Egg Stage Survival'})
subplot(3,2,6); contourf (EM, LM,ELAdultBeesW3); title({'AdultBees Number at Third End of Winter '});xlabel({'Larval Stage Survival'});ylabel({'Egg Stage Survival'})

set(h,'Position',[0.01 0.03 0.05 0.7])

contourf (EM, LM, ELM0x2DAdultBees1, 'DisplayName', 'EM, LM, ELM0x2DAdultBees1'); figure(gcf)


subplot(2,2,1); mesh(ELM0x2DAdultBees1)
subplot(2,2,2); mesh(ELM0x2DAdultBees2)



subplot(2,2,3); mesh(ELM0x2DHoneystore1)
subplot(2,2,4); mesh(ELM0x2DAdultBees1,ELM0x2DAdultBees2,ELM0x2DHoneystore1)

