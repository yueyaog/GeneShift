#!/usr/bin/env python

# Time Series Gene Expression Levels
# Gao's data. Working with univariate time series data.

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score
from tensorflow import keras
from sklearn.preprocessing import OneHotEncoder
from sklearn.metrics import confusion_matrix
import argparse


# Setup command line arguments
parser = argparse.ArgumentParser(description='classifier')
parser.add_argument('--file_name', type=str, default='my_file', help="The name of the folder ")
#parser.add_argument('--run_number', type=str, default='run1', help="The trial number of the run ")
parser.add_argument('--num_features', type=int, default=1, help="The number of features in the data set. Should only be 1 (time)")
parser.add_argument('--threshold', type=int, default=2, help="The minimum number of samples for each cluster ")
parser.add_argument('--test_split', type=float, default=0.3, help="The percentage of samples to split")
parser.add_argument('--epochs', type=int, default=100, help="Number of passes through the dataset")

args = parser.parse_args()

# Load in data
df = pd.read_csv(args.file_name)

# Create output folder to store results
OUTPUT_FOLDER = os.path.join("results",args.file_name)
if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)

# Create a log file 
text_file = open(OUTPUT_FOLDER + "/" + "output.txt", "a+")

# Remove clusters if they're not below a threshold
# Threshold must be >= 2 in order to create the train and test samples
label_counts = df.iloc[:,-1].value_counts()
index_counts = label_counts.index.tolist()
clusters = df.iloc[:,-1]
for count, index in zip(label_counts, index_counts):
    if count < args.threshold:
        text_file.write("\nCluster: %s was removed from analysis" % (index))
        df  = df[clusters != index]

# Extract data and labels/convert to one-hot
X = df.iloc[:,:-1] # data and gene name
y = df.iloc[:,-1]
y_labels = y.copy()
data_df = X.iloc[:,1:] # only data
num_labels = len(set(y_labels))

# Encode clusters into one-hot labels
onehot_encoder = OneHotEncoder()
y_labels = y_labels.values.reshape(-1, 1)
y_one_hot = onehot_encoder.fit_transform(y_labels)

# Grab the encoded classes
classes = onehot_encoder.categories_[0]

# Count of each label 
fig, ax = plt.subplots(figsize=(20,7.5))
#sns.countplot(x='Cluster', data=df)
sns.countplot(x=clusters, data=df)

ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
plt.tight_layout()
plt.savefig(OUTPUT_FOLDER + "/countplot.png", dpi=300)
plt.show()

# Generate random line graphs from the data
rows = 3
cols = 3
plot_size = 3.5

indices = np.random.choice(np.arange(len(data_df)), rows * cols)

plt.figure(figsize=(plot_size * cols, plot_size * rows))

for i in range(rows * cols):
    index = indices[i]
    
    ax = plt.subplot(rows, cols, i + 1)
    X.iloc[index,1:].plot(title=X.iloc[index][0])
    #plt.imshow(255 - X_train[index], cmap="gray")
    #plt.title("label = %d" % y[index])
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)

# Data is in the form of [samples, timesteps]
num_features_lstm = args.num_features
time_steps = data_df.shape[1]

# For LSTM need to reshape data to the form of [samples, timesteps, features]
data_df_lstm = data_df.values.reshape((data_df.shape[0], time_steps, num_features_lstm))
X_lstm = X.values.reshape((X.shape[0], X.shape[1], num_features_lstm))


# Split Data and Train Model
# The CNN and Linear model are 2D and can use the same split variables
# The LSTM model is now 3D data and requires different split variables

# Random state splits split the data the same for both data splits
X_train_lstm, X_test_lstm, y_train_lstm, y_test_lstm = train_test_split(data_df_lstm, y_one_hot, stratify=y, test_size = args.test_split, random_state = 42)
# Split linear model data
X_train, X_test, y_train, y_test = train_test_split(data_df, y_one_hot, stratify=y, test_size = args.test_split, random_state = 42)

X_train_cnn = X_train.values.reshape((X_train.shape[0], X_train.shape[1], 1))
X_test_cnn = X_test.values.reshape((X_test.shape[0], X_test.shape[1], 1))

# Model parameters 
epochs = args.epochs
batch_size = 32
num_features = X_train.shape[1]
num_features_cnn = X_train_cnn.shape[2]

def create_lstm(num_labels, num_features, time_steps, layer_units=100, dropout=0.2):
    """Create LSTM model"""
    # Last LSTM should have reutrn_sequences=False for classification
    # define LSTM model
    model = keras.models.Sequential()
    # try changing return_sequences to false
    #model.add(keras.layers.LSTM(units=100, dropout=0.2, activation='relu', return_sequences=True, input_shape=(time_steps, num_features)))
    model.add(keras.layers.LSTM(units=100, dropout=dropout, activation='relu', return_sequences=True, input_shape=(time_steps, num_features)))
    #model.add(keras.layers.Flatten())
    model.add(keras.layers.LSTM(units=50, dropout=dropout, activation='relu', return_sequences=False))
    model.add(keras.layers.Dense(units=num_labels, activation='softmax'))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=["accuracy"])
    
    return model


def create_mlp(input_size, num_labels):
    """Create MLP model"""
    # create a 3-layer neural network
    model = keras.models.Sequential()
    model.add(keras.layers.Dense(units=1024, activation="relu", input_shape=(input_size,)))
    model.add(keras.layers.Dropout(0.2))
    model.add(keras.layers.Dense(units=512, activation="relu"))
    model.add(keras.layers.Dropout(0.2))
    model.add(keras.layers.Dense(units=256, activation="relu"))
    model.add(keras.layers.Dropout(0.2))
    model.add(keras.layers.Dense(units=num_labels, activation="softmax"))

    #model.add(keras.layers.Dense(units=num_labels, activation="softmax"))
    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=["accuracy"])
    
    return model


def create_cnn(time_steps, num_features, num_labels):
    """Create 1D CNNl model"""
    model = keras.models.Sequential()
    model.add(keras.layers.Conv1D(filters=128, kernel_size=3, strides=2, activation="relu", input_shape=(time_steps, num_features)))
    model.add(keras.layers.MaxPooling1D(pool_size=2))
    model.add(keras.layers.Dropout(0.2))
    model.add(keras.layers.Flatten())
    model.add(keras.layers.Dense(units=num_labels, activation='softmax'))

    model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=["accuracy"])
    
    return model

model_cnn = create_cnn(time_steps, num_features_cnn, num_labels)
model_lstm = create_lstm(num_labels, num_features_lstm, time_steps,  layer_units=100, dropout=0.2)
model_linear = create_mlp(num_features, num_labels)

history_cnn = model_cnn.fit(X_train_cnn, y_train, epochs=epochs, batch_size=batch_size, validation_split=0.1)
history_lstm = model_lstm.fit(X_train_lstm, y_train_lstm, epochs=epochs, batch_size=batch_size, validation_split=0.1)
history_linear = model_linear.fit(X_train, y_train, epochs=epochs, batch_size=batch_size, validation_split=0.1)

def evaluate_model(history, title_acc="Training Accuracy", title_loss="Training Loss"):
    # plot the training accuracy
    plt.plot(history.history["acc"])
    plt.plot(history.history["val_acc"])
    plt.title(title_acc)
    plt.ylabel("Accuracy")
    plt.xlabel("Epoch")
    plt.legend(["Training", "Validation"], loc="upper left")
    plt.show()
    
    # plot the training loss
    plt.plot(history.history["loss"])
    plt.plot(history.history["val_loss"])
    plt.title(title_loss)
    plt.ylabel("Loss")
    plt.xlabel("Epoch")
    plt.legend(["Training", "Validation"], loc="middle")
    plt.show()


# Evaluate model
evaluate_model(history_lstm, title_acc="Training Accuracy CNN", title_loss="Training Loss CNN")
evaluate_model(history_lstm, title_acc="Training Accuracy LSTM", title_loss="Training Loss LSTM")
evaluate_model(history_linear, title_acc="Training Accuracy Linear", title_loss="Training Loss Linear")

# Predict raw network values
y_raw_predict_cnn = model_cnn.predict(X_test_cnn)
y_raw_predict_lstm = model_lstm.predict(X_test_lstm)
y_raw_predict_linear = model_linear.predict(X_test)

# Getting index prediction for confusion matrix
y_pred_linear = np.argmax(y_raw_predict_linear, axis=1)
y_pred_cnn = np.argmax(y_raw_predict_cnn, axis=1)
y_pred_lstm = np.argmax(y_raw_predict_lstm, axis=1)
y_test_lstm = np.argmax(y_test_lstm.toarray(), axis=1)
y_test = np.argmax(y_test.toarray(), axis=1)

def plot_cm(y_test, y_pred, classes, cbar='False', title="Confusion Matrix", normalize=0, dpi=320):
    """ Compute confusion matrix for the ground truth and predicted labels.
        if normalize = 1 the confusion matrix will work as a %, not number o
    """ 
    fig, ax = plt.subplots(figsize=(30,30))         # Sample figsize in inches
    cnf_matrix = confusion_matrix(y_test, y_pred)
    
    if normalize:
        # Normalize confusion matrix
        cnf_matrix = ((cnf_matrix.astype('float') / cnf_matrix.sum(axis=1)[:, np.newaxis])*100).astype(int)
    
    # plot a heatmap of the confusion matrix
    sns.heatmap(cnf_matrix, annot=True, fmt="d", square=True, xticklabels=classes, yticklabels=classes)
    plt.ylabel("Expected")
    plt.xlabel("Measured")
    plt.title(title)
    plt.savefig(OUTPUT_FOLDER + "/" + title + '.png', dpi=dpi)
    plt.show()
    
plot_cm(y_test, y_pred_cnn, classes, cbar="True", title="cm_cnn_scaled", normalize=1)
plot_cm(y_test, y_pred_linear, classes, cbar="True", title="cm_mlp_scaled", normalize=1)
plot_cm(y_test_lstm, y_pred_lstm, classes, cbar="True", title="cm_lstm_scaled", normalize=1)

text_file.write("CNN f1 score: %0.2f\n" % (f1_score(y_test, y_pred_cnn, average="weighted")))
text_file.write("Linear f1 score: %0.2f\n" % (f1_score(y_test, y_pred_linear, average="weighted")))
text_file.write("LSTM f1 score: %0.2f\n" % (f1_score(y_test_lstm, y_pred_lstm, average="weighted")))

#text_file.write("\n%0.2f " % (f1_score(y_test, y_pred_cnn, average="weighted")))
#text_file.write("%0.2f " % (f1_score(y_test, y_pred_linear, average="weighted")))
#text_file.write("%0.2f \n" % (f1_score(y_test_lstm, y_pred_lstm, average="weighted")))


text_file.close()
