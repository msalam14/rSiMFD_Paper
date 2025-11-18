# Necessary libraries are imported here
import keras.optimizers
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#import PIL as pil
import skimage.io as skio
import sklearn as skl
import sklearn.metrics as sklm
from random import sample
import pathlib
import random
from skimage.color import rgb2gray
import keras.backend as K
import gc
from keras.models import Model
from keras.layers import Dense
from keras import optimizers

# We import the data set from tensorflow and build the model there
import tensorflow as tf
from tensorflow.keras import layers, models, metrics
from tensorflow.keras.preprocessing import image

def class_sum(pred_score,actual,thresh=0.5):
    predY=np.where(pred_score>=thresh,1,0)
    res=np.where(actual==predY,1,0)
    return([res.mean(),res[np.where(actual==1)].mean(),res[np.where(actual==0)].mean()])



## Meta data
path = "/Users/samsul/Project2WSI/peso_testset_png/"
meta_file= "peso_testset_mapping.csv"
info_dat = pd.read_csv(path+meta_file)
info_dat['Cancer'] = np.where(info_dat['type'] == 'cancer', 1, 0)


# data frame for the results
ntrain=25
cnn_res=pd.DataFrame(columns=['Accuracy','Sensitivity','Specificity','tAccuracy','tSensitivity','tSpecificity'])

# model fit params
npch=5
bsize=10

inp_shp1=224
inp_shp2=224
dep_shp=3

# model fitting
for l in range(ntrain):
    # training data preparation
    train_images = np.full((120, inp_shp1, inp_shp2,dep_shp), 0, dtype='float')
    trainY = np.full((120, 1), 0)
    test_images = np.full((40, inp_shp1, inp_shp2,dep_shp), 0, dtype='float')
    testY = np.full((40, 1), 0)

    random.seed(10+l)
    tst_sam = random.sample(np.ndarray.tolist(pd.unique(info_dat["slide"])), 30)

    train_meta = info_dat[info_dat['slide'].isin(tst_sam)]
    test_meta = info_dat[~info_dat['slide'].isin(tst_sam)]

    # subject in the train set

    imID = train_meta['id'].tolist()

    # training formation

    for k in range(train_meta.shape[0]):
        fname = str(imID[k]) + '.jpg'
        fpath = path + fname
        train_images[k, :, :] = image.img_to_array(image.load_img(fpath,target_size=(inp_shp1, inp_shp2,dep_shp)))/255

    # training class

    trainY = train_meta['Cancer'].to_numpy().reshape((120, 1))

    # subject in the test data
    timID = test_meta['id'].tolist()

    # test data formation

    for k in range(test_meta.shape[0]):
        fname = str(timID[k]) + '.jpg'
        fpath = path + fname
        test_images[k, :, :] = image.img_to_array(image.load_img(fpath,target_size=(inp_shp1, inp_shp2,dep_shp)))/255


    # test label

    testY = test_meta['Cancer'].to_numpy().reshape((40, 1))

    # CNN architecture
    tf.keras.backend.clear_session()
    # Definition of network architecture

    #model = models.Sequential()

    #base_model = tf.keras.applications.ResNet50(include_top=False,
    #                                      weights='imagenet',
    #                                      input_shape=(224, 224, 3),
    #                                      pooling='avg')

    base_model = tf.keras.applications.InceptionV3(include_top=False,
                                                 weights='imagenet',
                                                 input_shape=(inp_shp1, inp_shp2, dep_shp),
                                                 pooling='avg')


    x_top = base_model.output
    x_out = Dense(1, name='output', activation='sigmoid')(x_top)
    model = Model(base_model.input, x_out)

    for each_layer in base_model.layers:
        each_layer.trainable = True

    # model.summary()

    # model compilation
    opt = tf.keras.optimizers.legacy.Adam(learning_rate=0.001)


    model.compile(optimizer=opt,
                  loss='binary_crossentropy',
                  metrics=['accuracy'])


    history = model.fit(train_images, trainY, epochs=npch,
                            batch_size=bsize,
                            validation_split=0.20)

    prd_scr_tr=model.predict(train_images).flatten()
    prd_scr_ts = model.predict(test_images).flatten()

    smr_tr=class_sum(pred_score=prd_scr_tr, actual=trainY.flatten(), thresh=0.5)

    smr_ts=class_sum(pred_score=prd_scr_ts, actual=testY.flatten(), thresh=0.5)

    cnn_res.loc[l]=[smr_tr[0],smr_tr[1],smr_tr[2],smr_ts[0],smr_ts[1],smr_ts[2]]

    del(train_images)
    del(trainY)



    #[loss, accr, prec, recl] = model.evaluate(x=test_images, y=testY, batch_size=5)

    #cnn_res.loc[l] = np.array([loss, accr, prec, recl])

    del(test_images)
    del(testY)
    del(history)
    gc.collect()


# saving the results

cnn_res.to_csv(path+'inceptionV3_train.csv',sep=',')



